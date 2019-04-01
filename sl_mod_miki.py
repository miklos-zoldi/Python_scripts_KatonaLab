# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 13:45:58 2014

@author: dudok
"""

"""
Barna's sl file modified by Miki (Jul 11 2018)

- if running in Pycharm
    1) one can uncomment lines #matplotlib.use('Qt4Agg')
    2) in selector2D/3D functions, instead of having spyder_blocking_bug lines, one can uncomment im.show(Block=True)
- if running in Spyder,
    1) comment out #matplotlib.use('Qt4Agg') line; adjust the default GUI to use Qt (Preferences -> Graphics)
    2) image.show() blocking behavior is not working, image is closed before it is displayed
    to overcome this, check spder_blocking_bug lines
    Present in every function using mpl_connect and in class draw_ROI key_press_callback (event==q/x, when one would close the plot)
    
Other: 
- getsub function: modified in order to handle 1 pointlist (1 ROI) or more pointlists (more ROIs) as inputs. Input should be list of list(s), e.g. [[1,2,3]]
- selector2D function: 
    modified to have similar output as selector3D and to be able to handle multiple ROIs
    lassocoordset is saved as self.lasso for being able to display it later
    optionally can add zooming
- selector3D function: 
    lassocoordset is saved as self.lasso for being able to display it later
    optionally can add zooming
- showstorm function: 
    freehand selection (lasso) can be displayed if specified
- slices function: 
    freehand selection (lasso) can be displayed if specified
    Stormpoints to be added to the MIP are named correctly
        x_zy ---> x_xz                                           
        y_xz ---> zy
    Using median of z data range instead of zmin/zmax when determining MIPrange              
- drawROI class:
    def zoom and things related to it (in def key_press_callback event='r'; in def __init__: self.orig_xlim/ylim )
- coords class + colums included 
- showclasts(): handling enumerate object
- class image: 
    handling imageDesc = tif[0].description error (TypeError:'TiffFile' object does not support indexing)


SHOULD CORRECT:
- closing the figure in selector2D/3D does not work, it is stucked in an infinite loop, next image not displayed
PRESS X INSTEAD if using Spyder
- selector3D: why thisfilt is not erased in every loop???
- making spyder vs pycharm user-friendly (Program would choose which lines should be considered)
- have to press keypress_callback keys twice
- getsub - error message if no point are included and 'q' was pressed
"""
import cfg
import numpy
import scipy
import math
from scipy import ndimage, signal, stats
from scipy.spatial import ConvexHull
import PIL.Image
import PIL.ImageDraw
import PIL.ImageFilter
import copy
import os
from time import time
import pylab
#import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.path as mplpath
import matplotlib.pyplot as plt
import matplotlib.cm as colormap
import cv2
import operator
from random import random, randint, randrange
try:
    import tifffile
    import xml.etree.ElementTree as XMLET
except:
    print 'tifffile or xml modules not found, methods with confocal images will not run'


class coords(object):
    def __init__(self, path=None, fname=None, Ch_ex='Z Rejected', Ch_in=None,
                 Fr=False, keeplines=False, thresh=False, accu=False,
                 rich=False, extracol=False):
        '''reads data if file is specified'''
        if rich:
            Fr=True
            keeplines=True
            accu=True
        if fname is None:
            self.nlp=0
            self.l=[]
            return
        xcoords=[]
        ycoords=[]
        zcoords=[]
        chs=[]
        frs=[]
        phots=[]
        acc=[]
        extra=[]
        ax=[]
        chnum=[]
        BG=[]
        frlengths=[]
        n_channels=0
        chids={}
        if keeplines:
            self.lines=[]
        self.xcol=22
        self.ycol=23
        self.zcol=17
        self.frcol=12
        self.chcol=0
        self.photonscol=18
        self.accucol=19
        self.axcol=9
        self.BGcol=10
        self.frlengthscol=13
        with open(path+fname) as f:
            lines=f.readlines()
        self.header=lines[0]
        for i in xrange(1, len(lines)):
            if lines[i]=='\n' or lines[i]=='\r\n':
                continue
            items=lines[i].split('\t')
            incl=True
            try:
                if not float(items[self.photonscol])>thresh:
                    incl=False
            except:
                pass
            try:
                if items[self.chcol] in Ch_ex:
                    incl=False
            except:
                pass
            try:
                if not items[self.chcol] in Ch_in:
                    incl=False
            except:
                pass
            if incl:
                xcoords.append(int(float(items[self.xcol])))
                ycoords.append(int(float(items[self.ycol])))
                zcoords.append(int(float(items[self.zcol])))
                c=items[self.chcol]
                if c in chids:
                    cn=chids[c]
                else:
                    chids[c]=n_channels
                    cn=n_channels
                    n_channels+=1
                chs.append(c)
                chnum.append(cn)
                if Fr:
                    frs.append(int(items[self.frcol]))
                if keeplines:
                    self.lines.append(lines[i])
                if thresh is True:
                    phots.append(int(float(items[self.photonscol])))
                if accu:
                    acc.append(float(items[self.accucol]))
                if extracol:
                    extra.append(items[extracol])
                if rich:
                    ax.append(float(items[self.axcol]))
                    BG.append(float(items[self.BGcol]))
                    frlengths.append(float(items[self.frlengthscol]))
        self.data=numpy.array((xcoords, ycoords, zcoords))
        self.channels=chs
        self.chnames=set(self.channels)
        self.chdict=chids
        self.chnum=chnum
        if Fr:
            self.frames=frs
        if thresh is True:
            self.phots=phots
        if accu:
            self.acc=acc
        if extracol:
            self.extracol=extracol
            self.extra=extra
        if rich:
            self.ax=ax
            self.BG=BG
            self.frlengths=frlengths
        self.tag=fname
        self.init2()
        
    def save(self, path, fname, incl=None, mockline=None, defheader=None):
        '''saves the storm data (with modified coordinates) of the optional list, or all molecules'''
        if mockline is None:        
            mockline='mock\t21675.9\t15206.8\t21675.9\t15206.8\t1683.32849\t13224.08496\t359.73749\t1.11619\t0.96793\t7193.72705\t14167.12500\t172\t1\t6\t1\t1937.17\t-249\t2380.33529\t15.78044\t21675.9\t15206.8\t21675.0\t15206.0'.split('\t')
        if defheader is None:
            defheader='Channel Name	X	Y	Xc	Yc	Height	Area	Width	Phi	Ax	BG	I	Frame	Length	Link	Valid	Z	Zc	Photons	Lateral Localization Accuracy	Xw	Yw	Xwc	Ywc'
        if incl is None:
            incl=self.l
        f=open(path+fname, 'w')
        try:
            f.write(self.header)
        except:
            f.write(defheader+'\n')
        try:
            xcol=self.xcol
            ycol=self.ycol
            zcol=self.zcol
            chcol=self.chcol
        except:
            chcol=0
            xcol=22
            ycol=23
            zcol=17
        for i in incl:
            try:
                items=self.lines[i].split('\t')
            except:
                #Use ONLY chcol, xcol, ycol, zcol. All the other values are get from the mockline and have a constant value          
                items=mockline
            items[chcol]=str(self.channels[i])
            items[xcol]=str(self.xc(i))
            items[ycol]=str(self.yc(i))
            items[zcol]=str(self.zc(i))
            textout='\t'.join(items)
            if textout[-1]!='\n':
                textout+='\n'
            f.write(textout)
        f.close()

    
    def fill(self, data, tag='noname', Ch=False, ChNum=False, Fr=False):
        '''fill empty object with the data in arguments'''
        self.data=data
        self.chdict={}
        if not Ch is False:
            self.channels=Ch
        else:
            self.channels=['Fill']*len(data[0])
        if not ChNum is False:
            self.chnum=ChNum
            for i, c in enumerate(set(ChNum)):
                self.chdict[c]=1
        else:
            self.chnum=[0]*len(data[0])
            self.chdict[0]=0
        if not Fr is False:
            self.frames=Fr
        else:
            self.frames=[0]*len(data[0])
        self.tag=tag
        self.init2()
        
    def add_channel(self, storm): 
        '''adds the storm data in the coords object to the current instance as a new channel'''
        self.data=numpy.append(self.data.transpose(),storm.data.transpose()).reshape((self.nlp+storm.nlp,3)).transpose()
        newch=len(self.chdict)
        self.chnum.extend([newch]*len(storm.data[0]))
        self.chdict[newch]=newch
        self.init2()
        for attr in ['acc','lines','frames','phots','channels','extra','ax']:
            if hasattr(self,attr) and hasattr(storm,attr):
                getattr(self,attr).extend(getattr(storm,attr))

    def oldxlsread(self, path, fname):
        '''reads old style imagej coordinates output file'''
        xcoords=[]
        ycoords=[]
        zcoords=[]
        chs=[]
        self.lines=[]
        self.xcol=3
        self.ycol=4
        self.zcol=5
        self.chcol=2
        roicol=6
        with open(path+fname) as f:
            lines=f.readlines()
        self.header=lines[0]
        for i in xrange(1, len(lines)):
            items=lines[i].split('\t')
            incl=True
            if items[roicol] == '0':
                incl=False
            if incl:
                xcoords.append(int(float(items[self.xcol])))
                ycoords.append(int(float(items[self.ycol])))
                zcoords.append(int(float(items[self.zcol])))
                chs.append(items[self.chcol])
        self.data=numpy.array((xcoords, ycoords, zcoords))
        self.channels=chs
        self.chnames=set(self.channels)
        self.tag=fname
        self.init2()
    
    def init2(self):
        self.nlp=len(self.data[0])
        if self.nlp>0:
            self.l=range(self.nlp)
            self.xmin=numpy.min(self.data[0])
            self.xmax=numpy.max(self.data[0])
            self.ymin=numpy.min(self.data[1])
            self.ymax=numpy.max(self.data[1])
            self.zmin=numpy.min(self.data[2])
            self.zmax=numpy.max(self.data[2])
            self.cmass=[numpy.mean(self.data[0]),numpy.mean(self.data[1]),
                        numpy.mean(self.data[2])]
    def xc(self, i):
        '''x coordinate of point i'''
        return self.data[0,i]
    def yc(self, i):
        '''y coordinate of point i'''
        return self.data[1,i]
    def zc(self, i):
        '''z coordinate of point i'''
        return self.data[2,i]
    def getc(self, i):
        '''list of coordinates of point i'''
        return [self.xc(i), self.yc(i), self.zc(i)]
    def dist(self,i,j):
        '''distance of points i and j'''
        x=self.xc(i)-self.xc(j)
        y=self.yc(i)-self.yc(j)
        z=self.zc(i)-self.zc(j)
        return int(math.sqrt(x*x+y*y+z*z))
        
    def center_z(self):
        '''modifies the z coordinates to have a mean of 0'''
        self.data[2,:]-=numpy.average(self.data[2,:])
    
    def offset_z(self, zoffset):
        '''adds the offset to the z coordinates'''
        self.data[2,:]+=zoffset
        
    def zhist(self, binsize=50):
        '''Display histogram of z coordinates. 
        optional: bin size'''
        x=self.data[2,:]
        nbins=int((max(x)-min(x))/binsize)
        plt.hist(x,bins=nbins)
        return plt
    
    def calc_dm(self):
        '''compute distance matrix and return all distances'''
        dm=numpy.empty((self.nlp,self.nlp), dtype=numpy.int)
        rkf=[]
        for i in self.l:
            for j in self.l[:i]:
                d=self.dist(i,j)
                dm[i,j]=d
                dm[j,i]=d
                rkf.append(d)
            dm[i,i]=0
        self.dm=dm
        return rkf
    
    def calc_3dhull(self, vol=False):
        '''compute 3D convex hull of all points'''
        self.hull3=scipy.spatial.ConvexHull(numpy.asarray(self.rs()), incremental=False, qhull_options=None)
        if vol:
            def tetrahedron_volume(a, b, c, d):
                return numpy.abs(numpy.einsum('ij,ij->i', a-d, numpy.cross(b-d, c-d))) / 6
            dt = scipy.spatial.Delaunay(self.rs())
            tets = dt.points[dt.simplices]
            vol = numpy.sum(tetrahedron_volume(tets[:, 0], tets[:, 1], 
                                tets[:, 2], tets[:, 3]))
            return vol
        
    def calc_2dhull(self):
        '''compute 3D convex hull of all points'''
        self.hull2=scipy.spatial.ConvexHull(numpy.asarray(self.rs()[:,:2]))
    
    def dist_pfromhull(self, point):
        '''measure distance of point from the convex hull.
        Point is expected as [x,y,z]'''
        
        try:
            self.hull3
        except:
            self.calc_3dhull()
        hull=self.hull3
        # Translate all points with vector 'point'
        trans_points = numpy.copy(hull.points)
        for k in range(len(trans_points)):
            trans_points[k][0] = trans_points[k][0] - point[0]
            trans_points[k][1] = trans_points[k][1] - point[1]
            trans_points[k][2] = trans_points[k][2] - point[2]

            # CALCULATING THE CLOSEST HULL DISTANCES, AND HULL SURFACE POINT
        distance = 100000000
        for triang in hull.simplices:  # Calculate the distance from all the hull triangles hull_points contains
        # the hull points
            constraint = 0  # 0 :If constrains satisfied 1: if constr not satisfied
            cround = 1
            x0 = [0.3, 0.3, 0.4]
            while constraint == 0:
                P1 = trans_points[triang[0]]
                P2 = trans_points[triang[1]]
                P3 = trans_points[triang[2]]

                def func(x):
                    """ Objective function """
                    return ((x[0] * P1[0] + x[1] * P2[0] + x[2] * P3[0]) ** 2 + (
                        x[0] * P1[1] + x[1] * P2[1] + x[2] * P3[1]) ** 2 + (
                                x[0] * P1[2] + x[1] * P2[2] + x[2] * P3[2]) ** 2 )

                cons = ({'type': 'eq', 'fun': lambda x: numpy.array([x[0] + x[1] + x[2] - 1])},
                        {'type': 'ineq', 'fun': lambda x: numpy.array([x[0]])},
                        {'type': 'ineq', 'fun': lambda x: numpy.array([x[1]])},
                        {'type': 'ineq', 'fun': lambda x: numpy.array([x[2]])})

                res = scipy.optimize.minimize(func, x0, constraints=cons, method='SLSQP', options={'disp': False})
                qtmp = res.x  # The optimized parameter values
                constraint = 1
                if qtmp[0] > 1.001 or qtmp[0] < -0.001 or qtmp[1] > 1.001 or qtmp[1] < -0.001 or qtmp[2] > 1.001 or \
                                qtmp[2] < -0.001:
                    constraint = 0
                    # Make a starting guess
                    x0a = randint(0, 90)
                    x0b = randint(0, 95 - x0a)
                    x0c = 100 - x0a - x0b
                    x0 = [float(x0a) / 100, float(x0b) / 100, float(x0c) / 100]

                # Select the smallest triangle distance
                Q = [0, 0, 0]
                Q[0] = qtmp[0] * P1[0] + qtmp[1] * P2[0] + qtmp[2] * P3[
                    0]  # The coordinates of the closest point on the hull (Q)
                Q[1] = qtmp[0] * P1[1] + qtmp[1] * P2[1] + qtmp[2] * P3[1]
                Q[2] = qtmp[0] * P1[2] + qtmp[1] * P2[2] + qtmp[2] * P3[2]
                distmp = math.sqrt(Q[0] ** 2 + Q[1] ** 2 + Q[2] ** 2)
                if distmp < distance and constraint == 1:
                    distance = distmp
                    qpoint_tr = Q
                    tri = triang
                cround = cround + 1
        qpoint = [qpoint_tr[0] + point[0], qpoint_tr[1] + point[1],
                  qpoint_tr[2] + point[2]]  # Closest point on the surface in the original coordinates
        return ([qpoint, distance, tri])   

        
    def dbscan(self, minN, eps):
        '''Performs DBScan clustering on the coordinates'''
        c=[]
        members=[]   
        noise=[]
        if self.nlp>minN:
            self.calc_tree()
            visited=[False]*self.nlp
            member=[False]*self.nlp
            for i in self.l:
                if not visited[i]:
                    visited[i]=True
                    friends=self.tree.query_ball_point(self.getc(i), eps)
                    if len(friends)>=minN:
                        nc=[]
                        c.append(nc)
                        nc.append(i)
                        members.append(i)
                        member[i]=True
                        for j in friends: 
                            if not visited[j]:      
                                visited[j]=True
                                jfriends=self.tree.query_ball_point(self.getc(j), eps)
                                if len(jfriends)>=minN:
                                    for k in jfriends:
                                        if not k in friends:
                                            friends.append(k)
                            if not member[j]:
                                nc.append(j)
                                members.append(j)
                                member[j]=True

            for i in self.l:
                if not member[i]:
                    noise.append(i)
        elif self.nlp>0:
            noise=self.l
        self.clusters=c
        self.cmembers=members
        self.cnoise=noise    
        if len(c)>0:
            self.meanclnum=float(len(members))/len(c)
        else:
            self.meanclnum=''
        return members
            
    def showclusts(self, noise=True):
        'Plots the clusters in 2D scatterplot. Including noise is optional.'''
        cl=self.clusters
        x=[]
        y=[]
        c=[]
        s=[]
        if len(cl)>0:
            colors=range(len(cl))
            for i, cl in enumerate(self.clusters):
                print i, cl
                ci=randrange(0,len(colors))
                color=colors[ci]
                colors.pop(ci)
                
                for j in cl:
                    x.append(self.xc(j))
                    y.append(self.yc(j))                   
                    c.append(color)                   
                    s.append(30)
        else:
            i=0  
        if noise:
            for j in self.cnoise:
                x.append(self.xc(j))
                y.append(self.yc(j))
                c.append(len(self.clusters)+1)
                s.append(2) 
        plt.scatter(x,y,s=s,c=c,alpha=0.5,cmap=colormap.hsv_r)
        plt.axes().set_aspect('equal')
        return plt
    
    def clsave(self, path, fname, noise=False):
        '''Saves the clusters in a file. Including noise is otional'''
        of=f(path,fname,'ClID,X,Y,Z')
        for i, cl in enumerate(self.clusters):
            for j in cl:
                of.w([i,self.xc(j),self.yc(j),self.zc(j)])
        if noise:
            for j in self.cnoise:
                of.w(['noise',self.xc(j),self.yc(j),self.zc(j)])
        of.cl()
                
    def mean_nn(self, eps=160):
        '''Returns the number of neighbors within eps, and stores the number of neighbors for individual localizations'''
        self.calc_tree()
        nn=[]
        for i in self.l:
            nn.append(self.tree_nn(i, eps))
        self.nn=nn
        return numpy.mean(nn)
        
    def lfunc(self,dlist,area):
        '''Returns the list of values of Ripley's L function for a list of distances'''
        self.calc_tree()
        d=max(dlist)
        sdm=numpy.zeros((self.nlp,self.nlp), dtype=numpy.int)
        for i in self.l:
            for j in self.tree.query_ball_point(self.getc(i), d):
                if j<i:
                    sdm[i,j]=self.dist(i,j)
        p=numpy.nonzero(sdm)
        f=math.pi*self.nlp*(self.nlp-1)*(self.nlp-1)
        lfv=[]
        for d in dlist:
            count=1
            for x in xrange(len(p[0])):
                if sdm[p[0][x],p[1][x]]<d:
                    count+=1
            lfv.append([d,math.sqrt((area*count)/f)])
        return lfv
        
    def lfuncn(self,dlist,area):
        '''Returns the list of neighbor counts of Ripley's L function for a list of distances'''
        self.calc_tree()
        d=max(dlist)
        sdm=numpy.zeros((self.nlp,self.nlp), dtype=numpy.int)
        for i in self.l:
            for j in self.tree.query_ball_point(self.getc(i), d):
                if j<i:
                    sdm[i,j]=self.dist(i,j)
        p=numpy.nonzero(sdm)
        lfv=[]
        for d in dlist:
            count=0
            for x in xrange(len(p[0])):
                if sdm[p[0][x],p[1][x]]<d:
                    count+=1
            lfv.append([d,count])
        return lfv
        
    def EucDistance2D(self, point1, point2):
        # simple Euclidean distance calculation in 2D
        dist = math.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)
        return dist
    
    def TriangleArea2D(self, point1, point2, point3):
        '''calculates area of the triangle of 2D points'''
        a = self.EucDistance2D(point1, point2)
        b = self.EucDistance2D(point2, point3)
        c = self.EucDistance2D(point3, point1)
        s = (a + b + c) / 2
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        return area
        
    def hull2darea(self):
        '''Returns the area of the 2D convex hull of the points'''
        hull2D=ConvexHull(self.rs2d())
        hull2D_area = 0.0
        hull2D_peri = 0
        xpoints=[]
        ypoints=[]
        for k in range(len(hull2D.vertices) - 2):
            
            hull2D_area += self.TriangleArea2D(hull2D.points[hull2D.vertices[-1]], hull2D.points[hull2D.vertices[k]],
                                               hull2D.points[hull2D.vertices[k + 1]])
        for k in range(len(hull2D.vertices)):
           xpoints.append(hull2D.points[hull2D.vertices[k]][0])
           ypoints.append(hull2D.points[hull2D.vertices[k]][1])                                 
        self.hpoints=[xpoints,ypoints]
        p=hull2D.points
        for edge in hull2D.simplices:
            hull2D_peri+=math.sqrt(
            (p[edge[0]][0]-p[edge[1]][0])**2+
            (p[edge[0]][1]-p[edge[1]][1])**2)
        self.hullperi=hull2D_peri
        return hull2D_area
                           
    def measure_inter(self):
        '''returns mean d/r ratio, stores individual values'''
        hull2D=ConvexHull(self.rs2d())
        hull2D_area = 0.0
        for k in range(len(hull2D.vertices) - 2):
            hull2D_area += self.TriangleArea2D(hull2D.points[hull2D.vertices[-1]], hull2D.points[hull2D.vertices[k]],
                                               hull2D.points[hull2D.vertices[k + 1]])
        r = math.sqrt(hull2D_area / math.pi)
        c=[numpy.mean(self.data[0]),numpy.mean(self.data[1]),numpy.mean(self.data[2])]
        ratio=numpy.empty(self.nlp)
        for i in self.l:
            x=self.xc(i)-c[0]
            y=self.yc(i)-c[1]
            z=self.zc(i)-c[2]
            ratio[i]=math.sqrt(x*x+y*y+z*z)/r
        self.ratio=ratio
        return numpy.mean(ratio)
    
    def tree_nn(self, i, eps):
        '''number of neigbors of point i within eps from kD tree'''
        return len(self.tree.query_ball_point(self.getc(i), eps))-1
        
    def tree_nndist(self, i, mineps, maxeps):
        '''Distance of nearest neighbor of point i from kD tree'''
        mindist=maxeps
        for j in self.tree.query_ball_point(self.getc(i), maxeps):
            if j!=i:
                d=self.dist(i,j)
                if d>mineps:
                    mindist=min(d, mindist)
        return mindist
        
    def nndist(self, mineps, maxeps):
        '''returns all nearest neighbor distances'''
        d=[]
        try:
            self.tree
        except:
            self.calc_tree()
        for i in self.l:
            d.append(self.tree_nndist(i, mineps, maxeps))
        return d
        
    def tree2d_nn(self, x, y, eps):
        '''number of points within eps from x,y'''
        return len(self.tree2d.query_ball_point([x,y], eps))
    
    def rs(self):
        '''returns restructured data as list of x-y-z triplets'''
        rsd=numpy.empty((self.nlp, 3), dtype=numpy.int)
        for i in self.l:
            rsd[i]=self.getc(i)
        return rsd
        
    def rs2d(self):
        '''returns restructured data as list of x-y pairs'''
        rsd=numpy.empty((self.nlp, 2), dtype=numpy.int)
        for i in self.l:
            rsd[i]=[self.xc(i),self.yc(i)]
        return rsd
        
    def calc_tree(self):
        '''compute 3D tree'''
        self.tree=scipy.spatial.cKDTree(self.rs())
        
    def calc_2dtree(self):
        '''compute 2D tree'''
        self.tree2d=scipy.spatial.cKDTree(self.rs2d())
        return self.tree2d
        
    def dens_filt(self, n, eps):
        '''filter data for localizations with n neighbors within eps'''
        filt=[]
        try:
            self.tree
        except:
            self.calc_tree()
        for i in self.l:
            if not self.tree_nn(i, eps)<n:
                filt.append(i)        
        self.filt=filt
        self.fnlp=len(filt)
        return filt
        
    def trim_z(self, limit):
        '''Returns list of molecules within given range from z=0'''
        filt=[]
        for i in self.l:
            if abs(self.zc(i))<limit:
                filt.append(i)
        self.filt=filt
        self.fnlp=len(filt)
        return filt
        
    def pick_params(self, wnlp):
        '''increases filtering strictness until the nlp is below the specified limit'''
        nset=(2,3,4)
        epset=(150,125,100,80,60,50)
        for n in nset:
            for eps in epset:
                self.dens_filt(n, eps)
                if self.fnlp<wnlp:
                    print '('+str(n)+', '+str(eps)+')',
                    return n, eps
        print '('+str(n)+', '+str(eps)+')',
        return n, eps
    
    def getsub(self, l):
        '''returns new instance containing subset of coordinates included in l'''
        data=numpy.empty((3, len(l)))
        for i,j in enumerate(l):
            data[:,i]=self.data[:,j]           
        nd=coords()
        nd.fill(data, tag=self.tag+' subset')
        for attr in ['acc','lines','frames','phots','channels','extra','ax','chnum','BG','frlengths']:
            if hasattr(self,attr):
                setattr(nd,attr,[getattr(self,attr)[i] for i in l])
        if hasattr(nd,'channels'):
            nd.chnames=set(nd.channels)
        nd.chdict=self.chdict
        return nd
        
    def getchannel(self,c):
        '''returns the selected channel specified by
        either channel name or channel number'''
        if type(c)==str:
            l=[]
            for i in self.l:
                if self.channels[i] == c:
                    l.append(i)
            return self.getsub(l)
        elif type(c)==int:
            l=[]
            for i in self.l:
                if self.chnum[i] == c:
                    l.append(i)
            return self.getsub(l)
        else:
            print 'str for channel name or int for channel number is expected'
    
    def gamma(self, image, g=2.2):
        '''apply gamma on image'''
        image[:,:]=((image[:,:]/255)**(1.0/g))*255
        return image

    def filt_fr(self, startf, stopf):
        '''Returns list of molecules included in specified frame range'''
        filt=[]
        for i in self.l:
            if self.frames[i]>=startf and self.frames[i]<stopf:
                filt.append(i)
        self.filt=filt
        return filt
    
    def count_frame(self, fr_start, fr_range):
        n=0
        for i in self.l:
            if self.frames[i]>=fr_start:
                n+=1
                if self.frames[i]>=fr_start+fr_range:
                    break
        return n
        
    def pix_show(self, psize=80, norm=True, dfilt=False, show=False, rdat=False, fullsize=False):
        '''display image of storm data'''
        if fullsize is True:
            fullsize=45000
        if fullsize:
            sx=int(fullsize/psize)
            sy=int(fullsize/psize)
            xmin=0
            ymin=0
        else:
            sx=int(self.xmax/psize)+1-int(self.xmin/psize)
            sy=int(self.ymax/psize)+1-int(self.ymin/psize)
            xmin=self.xmin
            ymin=self.ymin
        pim=numpy.zeros((sx,sy), dtype=numpy.float)
        if dfilt:
            l=self.filt
        else:
            l=self.l
        for i in l:
            x=int((self.xc(i)-xmin)/psize)
            y=int((self.yc(i)-ymin)/psize)
            pim[x,y]+=norm
        if norm is True:
            mv=numpy.amax(pim)
            pim*=255.0/mv
            pim=self.gamma(pim)
        self.pixim=PIL.Image.fromarray(numpy.array(pim, dtype=numpy.uint8))
        if show:
            self.pixim.show()
        if rdat:
            return pim
    
    def translate(self, xcorr, ycorr, apply_shift=False):
        '''modifies the x-y coordinates, or returns a modified copy'''
        if not apply_shift:
            data=copy.deepcopy(self.data)
            for i in self.l:
                data[0,i]+=xcorr
                data[1,i]+=ycorr
            nd=coords()
            nd.fill(data)
            return nd
        for i in self.l:
            self.data[0,i]+=xcorr
            self.data[1,i]+=ycorr
        
    def hotspot(self, r):
        '''identifies the center of a circle with maximal labeling intensity'''
        try:
            self.tree2d
        except:
            self.calc_2dtree()
        maxn=0
        step=int(r*0.25)
        for x in xrange(int(self.xmin), int(self.xmax)+step+1, step):
            for y in xrange(int(self.ymin), int(self.ymax)+step+1, step):
                n=self.tree2d_nn(x, y, r)
                if n>maxn:
                    maxn=n
                    maxpos=[x,y]
        self.hspos=maxpos
        return self.tree2d.query_ball_point(maxpos, r)
         
                
    def display(self, psize=80):
        '''shows image'''
        self.pix_show(psize, show=True)
                
    def __getitem__(self, key):
        return self.getc(key)
    def __iter__(self):
        return iter(self.rs())        
    def __str__(self):
        try:
            return 'STORM data from '+self.tag+':'+str(self.nlp)+' NLP'
        except:
            return 'empty'


class f(object):
    def __init__(self, basepath, name, vs, dumpy=False):
        '''Creates tab separated output file. Provide column names separated by comma.
        If dumpy, force writes each line to disk, otherwise writes everything at closing.'''
        self.dumpy=dumpy
        if not '.' in name:
            name+='.txt'
        self.f=open(basepath+name, 'w')
        self.fields=vs.split(',')
        self.header=vs.replace(',','\t')
        self.f.write(self.header+'\n')
        self.s=''
        
    def w(self, vs):
        '''Adds a tab separated line from list of values'''
        self.s+='\t'.join(str(v) for v in vs)+'\n'
        if self.dumpy:
            self.dump()
        
    def dump(self):
        self.f.write(self.s)
        self.f.flush()
        os.fsync(self.f.fileno())
        self.s=''
    
    def cl(self):
        '''Write to disc and close'''
        self.f.write(self.s)
        self.f.close()
        
        
class image(object):
    '''reads specified tif image with metadata'''
    def __init__(self, path, fname):
        self.file_path=path+fname
        self.isParsingNeeded = False
        self.ConfocalMetaData = {}
        with tifffile.TiffFile(self.file_path) as tif:
            # get metadata
            try:
                imageDesc = tif.pages[0].image_description
            except:
                im = PIL.Image.open(self.file_path)
                imageDesc = im.tag[270]
                imageDesc = str(imageDesc).strip("',)")[3:]
            root = XMLET.fromstring(imageDesc)
            for L1 in root:
                if L1.tag.split('}')[1] == 'Image':
                    for L2 in L1:
                        if L2.tag.split('}')[1] == 'Pixels':
                            if float(L2.attrib['PhysicalSizeX']) !=1.0:
                                self.ConfocalMetaData['SizeC'] = float(L2.attrib['SizeC'])
                                self.ConfocalMetaData['SizeX'] = float(L2.attrib['PhysicalSizeX'])
                                self.ConfocalMetaData['SizeY'] = float(L2.attrib['PhysicalSizeY'])
                                self.ConfocalMetaData['SizeZ'] = float(L2.attrib['PhysicalSizeZ'])
            #read image data
            self.ConfocalData = tif.asarray()
        # confocaldata first dimensions z channel, color channels, than data
        ConfShape = self.ConfocalData.shape
        ConfDim = len(self.ConfocalData.shape)
        if ConfDim == 4:
            self.ConfocalMetaData['ChannelNum'] = int(ConfShape[1])
        if ConfDim == 3:
            if 'SizeC' in self.ConfocalMetaData:
                self.ConfocalMetaData['ChannelNum'] = int(self.ConfocalMetaData['SizeC'])
            else:
                self.ConfocalMetaData['ChannelNum'] = 1
        if ConfDim == 2:
            self.ConfocalMetaData['ChannelNum'] = 1
        #set the z channel to the middle of the stacks
        NumOfZSlices=0
        # more Zchannels and color channels
        if len(self.ConfocalData.shape) > 3:
                NumOfZSlices=self.ConfocalData.shape[0]           
        elif len(self.ConfocalData.shape) > 2 and self.ConfocalMetaData['ChannelNum'] == 1:
                #z channels only
                NumOfZSlices=self.ConfocalData.shape[0]       
        elif len(self.ConfocalData.shape) > 2 and self.ConfocalMetaData['ChannelNum'] > 1:
                #color channels only
                NumOfZSlices=1
        else:
                #a single tif image
                NumOfZSlices=1

class bouton(object):
    '''collects all data from multiple files about the analyzed bouton'''
    def __init__(self, wdir, tag, res='_Results.txt', roicoord='_RoiCoords.txt', roi='_RoiAttr.txt', cf=None):
        if '_Results.txt' in tag:
            self.tag=tag[:tag.find('_Results.txt')]
        else:
            self.tag=tag
        self.wdir=wdir
        self.cell=self.tag.split('_')
        if res:
            self.getres(res)
        if roicoord:
            self.getroistorm(roicoord)
        if roi:
            self.getroi(roi)
        if cf:
            self.getcf(cf)
                
    def getres(self, ext):
        '''reads vividstorm results file into a dict'''
        f=self.tag+ext
        res={}
        try:
            with open(self.wdir+f) as fi:
                lines=fi.readlines()
            k=lines[0][:-1].split('\t')
            v=lines[1][:-1].split('\t')
            for i in xrange(len(k)):
                nv=v[i]
                try:
                    nv=float(nv)
                except:
                    pass
                key=k[i]
                res[key]=nv
        except:
            print 'Results file not found:'+f
        self.res=res

    def getroistorm(self, ext):
        '''reads coordinates file into coords object'''
        try:
            self.storm=coords(self.wdir, self.tag+ext, rich=True)
        except:
            print 'Coordinates file not found:'+self.tag+ext
            self.storm=None
        
    def getroi(self, ext):
        '''reads vividstorm roi attributes file'''
        f=self.tag+ext
        try:     
            with open(self.wdir+f) as fi:
                lines=fi.readlines()                
            l=lines[1].split('\t')
            self.offset=[float(l[0]),float(l[1][:-1])]
            l=lines[7][:-1].split('\t')[:-1]
            self.cfint=[float(v) for v in l]
            self.cfname=lines[9][:-1]
            roix=[]
            roiy=[]
            for i in xrange(11,len(lines)):
                l=lines[i].split('\t')
                roix.append(float(l[0]))
                roiy.append(float(l[1]))
            self.roi=[roix,roiy]
        except:
            print 'ROI import failed:'+f
            self.offset=None
            self.roi=None
            
    def roipoly(self):
        '''returns a polygon object of the roi'''
        l=len(self.roi[0])
        r=numpy.empty((l,2))
        for i in xrange(l):
            r[i,:]=[self.roi[0][i],self.roi[1][i]]
        return mplpath.Path(r)
            
    def roiperi(self):
        '''returns perimeter of roi in microns'''
        x=self.roi[0]
        y=self.roi[1]
        p=0
        for i in xrange(1,len(x)):
            p+=math.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)
            if int(x[0])==int(x[i]) and int(y[0])==int(y[i]):
                break
        return p/1000
            
    def getcf(self,path):
        '''finds confocal image belonging to the bouton'''
        if not path.endswith('//'):
            path+='//'
        fn=self.cfname
        self.cf=image(path, fn)
        self.calib=self.cf.ConfocalMetaData['SizeX']*1000
    
    def get_threshold(self,c,w):
        '''returns threshold for channel c in slice z from the center of size w*w'''
        start=numpy.size(self.cf.ConfocalData[0,0,0,:])/2-w/2
        z=len(self.cf.ConfocalData)
        thrs=[]
        for zs in xrange(z):
            data=numpy.empty((w,w), dtype=numpy.uint8)        
            data[:,:]=self.cf.ConfocalData[zs,c,start:start+w,start:start+w]
            ret1=int(numpy.mean(data))+1
            ret,th = cv2.threshold(data,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
            thrs.append(ret+ret1)
        return max(thrs)
        
    def pick_random(self,n,zrange):
        '''returns a coords object containing n points at random
        positions inside the roi, within +-zrange'''
        xmin,xmax=min(self.roi[0]),max(self.roi[0])
        ymin,ymax=min(self.roi[1]),max(self.roi[1])
        w,h=xmax-xmin,ymax-ymin
        poly=self.roipoly()        
        need=n
        data=numpy.empty((3,need))
        while need>0:
            x,y=xmin+random()*w,ymin+random()*h
            if poly.contains_point([x,y]):
                z=random()*zrange*2-zrange
                need-=1
                data[:,need]=[x,y,z]
        rnd=coords()
        rnd.fill(data)
        return rnd
        
    def fitroi(self, res):
        '''calculates all points of the roi with spline fit
        with set density of points (nm)'''
        x,y=self.roi[0],self.roi[1]
        origl=len(x)
        x.extend(x[:5])
        y.extend(y[:5])
        roi=numpy.array([x,y])
        x,y=roi[0],roi[1]
        t = numpy.zeros(x.shape)
        t[1:] = numpy.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2)
        t = numpy.cumsum(t)
        t /= t[-1]
        npoints=int(self.roiperi()*1000/res)+1
        nt = numpy.linspace(0, 1, npoints)
        x2 = scipy.interpolate.spline(t, x, nt)
        y2 = scipy.interpolate.spline(t, y, nt)
#        plt.plot(x2, y2, label='dist_spline')
#        assert False
        return [x2,y2]
        
    def showroi(self,channels=[[0,0],[1,1],[2,2]]):
        '''creates image of confocal image belonging to the area defined by vividstorm roi.
        channels parameter defines up to 3 channels: list of list from ch to ch. Example:
        [[0,0],[1,1],[2,2]] sends first 3 channels to r,g,b.'''
        xmin=(min(self.roi[0])-self.offset[0])/self.calib-3
        xmax=(max(self.roi[0])-self.offset[0])/self.calib+3
        ymin=(min(self.roi[1])-self.offset[1])/self.calib-3
        ymax=(max(self.roi[1])-self.offset[1])/self.calib+3
        w=int(xmax)-int(xmin)+1
        h=int(ymax)-int(ymin)+1
        im=numpy.zeros((h,w,3), dtype=numpy.uint8)
#        for c in channels[:min(3,len(self.cf.ConfocalData))]:
#            im[:,:,c[0]]=numpy.asarray(self.cf.ConfocalData[c[1]][int(ymin):int(ymax)+1,int(xmin):int(xmax)+1])
        self.im=PIL.Image.fromarray(im, 'RGB')
        self.bounding=[xmin,xmax+1,ymin,ymax+1]
        
    def showstorm(self,channels=[[0,0],[1,1],[2,2]],res=5,stormsize=10,
                  hideroi=False,hideconfoc=False,hidelasso=True,stormtype='outline',
                  colors=[(255,0,255),(0,255,255),(255,255,0)],
                  thr=None, blur=False, thstorm=None,fullimg=False,
                  resample=PIL.Image.BICUBIC):
        '''creates 2d view of the data'''
        stormsize=max(stormsize,res)
        s=self.storm
        xmin=(min(self.roi[0])-self.offset[0])/self.calib-3
        xmax=(max(self.roi[0])-self.offset[0])/self.calib+3
        ymin=(min(self.roi[1])-self.offset[1])/self.calib-3
        ymax=(max(self.roi[1])-self.offset[1])/self.calib+3
        w=int(xmax)-int(xmin)+1
        h=int(ymax)-int(ymin)+1
        if fullimg:
            xmin,ymin=128,128
            xmax,ymax=255+128,255+128
            h,w=256,256
        resc=self.calib/res
        self.showstorm_storm_xy=numpy.empty((s.nlp,2))
        im=numpy.zeros((h,w,3), dtype=numpy.uint8)
        if hideconfoc is False:
            hideconfoc=[False]*min(3,len(self.cf.ConfocalData))
        for i,c in enumerate(channels[:min(3,len(self.cf.ConfocalData))]):
            if hideconfoc[i]:
                continue
            im[:,:,c[0]]=numpy.asarray(
            self.cf.ConfocalData[c[1]][int(ymin):int(ymax)+1,
                                        int(xmin):int(xmax)+1])
        im=PIL.Image.fromarray(im,'RGB').resize((int(w*resc),int(h*resc)),resample=resample)
#        draw=PIL.ImageDraw.Draw(im)

        r=int(stormsize/res)
        if not thr is None and not thstorm:
            stormin,stormout=[],[]
        if blur:
            im=im.filter(PIL.ImageFilter.GaussianBlur(radius=blur))
        for i in s.l:
            x=int((float(s.xc(i)-self.offset[0])/self.calib-xmin)*resc+0.5)
            y=int((float(s.yc(i)-self.offset[1])/self.calib-ymin)*resc+0.5)
            self.showstorm_storm_xy[i,:]=[x,y]
            draw=PIL.ImageDraw.Draw(im)
            if stormtype=='outline':
                draw.ellipse((x-r, y-r, x+r, y+r), outline=colors[s.chnum[i]])
            if stormtype=='fill':
                draw.ellipse((x-r, y-r, x+r, y+r), fill=colors[s.chnum[i]])
            if not thr is None and not thstorm:
                if im.getpixel((x,y))[thr[0]]>thr[1]:
                    stormin.append(i)
                else:
                    stormout.append(i)
        if not thr is None and not thstorm:
            self.stormin=s.getsub(stormin)
            self.stormout=s.getsub(stormout)
        if not thstorm is None:      
            pdat=numpy.asarray(self.cf.ConfocalData[thr[0]][int(ymin):int(ymax)+1,
                                        int(xmin):int(xmax)+1])
            _,th=cv2.threshold(pdat,thr[1],1,cv2.THRESH_BINARY)
            nz=numpy.ndarray.nonzero(th)
            elem=numpy.transpose(nz)
            newp=numpy.empty((3,len(nz[1])))
            for i in xrange(len(nz[1])):
                y,x=[elem[i,0],elem[i,1]]
                xc=(x+xmin+0.5)*self.calib+self.offset[0]
                yc=(y+ymin+0.5)*self.calib+self.offset[1]                
                newp[:,i]=[xc,yc,0]
            thstorm=coords()
            thstorm.fill(newp)
            thstorm=self.roifilt(thstorm)
            self.thstorm=thstorm                   
        if not hideroi:
            pg=[] 
            for i in xrange(len(self.roi[0])):
                x=int((float(self.roi[0][i]-self.offset[0])/self.calib-xmin)*resc+0.5)
                y=int((float(self.roi[1][i]-self.offset[1])/self.calib-ymin)*resc+0.5)
                pg.append((x,y))
            draw.polygon(pg)
        if not hidelasso:
            lasso_numbers=len(self.lasso)
            for n in range(lasso_numbers):
                lasso_poly=[]
                for i in self.lasso[n]:
                    lasso_poly.append((i[0],i[1]))
                draw.polygon(lasso_poly)               
        self.im=im
        
    def crop(self,h,w,c,z,xmin,ymin):
        cropped=numpy.zeros((z,c,h,w), dtype=numpy.uint8)
        cropped[:,:c,:,:]=self.cf.ConfocalData[:,:c,ymin:ymin+h,xmin:xmin+w]
        return cropped
    
    def slices(self,zoffset=0, poly=False, res=5, resample=PIL.Image.BICUBIC,
               stormsize=10, hidestorm=False,hidelasso=True, stormdata=None, measure=None,
               thr=None, apply_thr=False, roifilt=True, hideroi=False,
               colors=[(255,0,255),(0,255,255),(255,255,255)],
               fullbounds=False,zmiprange=None,stormtype='outline'):
        '''create slices view of MIPs of the STORM data range of confocal ROI
        Expects full z stack'''
        s=stormdata
        if s is None:
            s=self.storm
        #crop roi
        xmin=(min(self.roi[0])-self.offset[0])/self.calib-3
        xmax=(max(self.roi[0])-self.offset[0])/self.calib+3
        ymin=(min(self.roi[1])-self.offset[1])/self.calib-3
        ymax=(max(self.roi[1])-self.offset[1])/self.calib+3
        w=int(xmax)-int(xmin)+1
        h=int(ymax)-int(ymin)+1
        z=len(self.cf.ConfocalData)
        c=min(3,self.cf.ConfocalMetaData['ChannelNum'])
        roi=self.crop(h,w,c,z,int(xmin),int(ymin))
        zcalib=self.cf.ConfocalMetaData['SizeZ']*1000
        resc=self.calib/res
        aspect=float(zcalib)/self.calib
        midslice=int((len(self.cf.ConfocalData)+1)/2)
        #save attribs for later reverse transform
        self.slices_start_zy=int((w+1)*resc+0.5)
        self.slices_start_xz=int((h+1)*resc+0.5)
        self.slices_midslice=midslice
        self.slices_xmin=xmin
        self.slices_ymin=ymin
        self.slices_resc=resc
        self.slices_aspect=aspect
        self.slices_zoffset=zoffset
        self.slices_zcalib=zcalib
        self.slices_w=w
        self.slices_h=h
        
        #threshold image
        if not thr is None:
            pdat=copy.deepcopy(roi)
            for ch,threshold in thr:
                arr=roi[:,ch,:,:]
#                blur=cv2.GaussianBlur(arr,(1,1),0)
                _,th=cv2.threshold(arr,threshold,1,cv2.THRESH_BINARY)
                scipy.ndimage.binary_erosion(th).astype(th.dtype)
                th*=255
                pdat[:,ch,:,:]=th[:,:,:]
                    
            #create threshold fake storm:
            if len(thr)>0:                
                points=numpy.empty((numpy.size(roi),3))
                j=0
                chs=[]
                for ch,_ in thr:
                    chs.append(ch)
                for ch in xrange(c):
                    if not ch in chs:
                        pdat[:,ch,:,:]=0
                nz=numpy.ndarray.nonzero(pdat[:,chs[0],:,:])
                elem=numpy.transpose(nz)
                for i in xrange(len(nz[1])):
                    points[j,:]=[elem[i,1],elem[i,2],elem[i,0]]
                    j+=1
                points=points[:j]
                for ch in chs[1:]:
                    newp=numpy.empty((j,3))
                    j=0
                    for x,y,zs in points:                      
                        if pdat[zs,ch,x,y]>0:
                            newp[j,:]=[x,y,zs]
                            j+=1
                    points=copy.deepcopy(newp[:j])
                thstorm=coords()
                if j>0:                    
                    newp=numpy.empty((3,j))
                    for j,p in enumerate(points):          
                        y,x,zs=p
                        xc=(x+self.slices_xmin+0.5)*self.calib+self.offset[0]
                        yc=(y+self.slices_ymin+0.5)*self.calib+self.offset[1]
                        zc=(zs-self.slices_midslice+0.5)*self.slices_zcalib-self.slices_zoffset
                        newp[:,j]=[xc,yc,zc]               
                    thstorm.fill(newp)
                    if roifilt:
                        thstorm=self.roifilt(thstorm)
                    if thstorm.nlp>0:
                        thstorm=thstorm.getsub(thstorm.dens_filt(4,200))
#                   s.add_channel(thstorm)
                    if thstorm.nlp>0:
                        s=thstorm
                self.thstorm=thstorm                
                if apply_thr:
                    roi=pdat
                
        #create mips
        self.roi_imgdata=roi 
        if s.nlp>3 and not fullbounds: 
            zstart=max(0,int(s.zmin/zcalib)+midslice-2)
            zstop=min(z,int(s.zmax/zcalib)+midslice+2)
#            zstart=max(0,int(numpy.median(s.data[2])/zcalib)+midslice-2)
#            zstop=min(z,int(numpy.median(s.data[2])/zcalib)+midslice+2)
            ystart=max(0,int((s.xmin-min(self.roi[0]))/self.calib)-2+3)
            ystop=min(w,int((s.xmax-min(self.roi[0]))/self.calib)+2+3)
            xstart=max(0,int((s.ymin-min(self.roi[1]))/self.calib)-2+3)
            xstop=min(h,int((s.ymax-min(self.roi[1]))/self.calib)+2+3)
            
        else:
            zstart,zstop=3,z-3
            ystart,ystop=3,h-3
            xstart,xstop=3,w-3
        if not zmiprange is None:
            zmipstart,zmipstop=max(0,midslice-zmiprange),min(z,midslice+zmiprange)
            zstart,zstop=zmipstart,zmipstop
            z=2*zmiprange+1            
        else:
            zmipstart,zmipstop=0,z
        
        xymip=numpy.zeros((h,w,c), dtype=numpy.uint8)
        for ch in xrange(c):
            for y in xrange(w):
                for x in xrange(h):
                    xymip[x,y,ch]=numpy.amax(roi[zstart:zstop+1,ch,x,y])
        
        xzmip=numpy.zeros((h,z,c), dtype=numpy.uint8)
        for ch in xrange(c):
            for zs in xrange(z):
                for x in xrange(h):
                    xzmip[x,zs,ch]=numpy.amax(roi[zs+zmipstart,ch,x,ystart:ystop+1])
                                       
        zymip=numpy.zeros((z,w,c), dtype=numpy.uint8)
        for ch in xrange(c):
            for y in xrange(w):
                for zs in xrange(z):
                    zymip[zs,y,ch]=numpy.amax(roi[zs+zmipstart,ch,xstart:xstop+1,y])
                         
        #assemble slices view
        slicesview=numpy.zeros((int(h*resc+z*resc*aspect+resc)+1,int(w*resc+z*resc*aspect+resc)+1,c), dtype=numpy.uint8)
        slicesview[:,:,:]+=255
        im=PIL.Image.fromarray(slicesview, 'RGB')
        im.paste(PIL.Image.fromarray(xymip,'RGB').resize((int(w*resc),int(h*resc)),
                         resample=resample),
                         box=(0,0))
        im.paste(PIL.Image.fromarray(xzmip,'RGB').resize((int(z*resc*aspect),int(h*resc)),
                         resample=resample),
                         box=(int(w*resc+resc),0))
        im.paste(PIL.Image.fromarray(zymip,'RGB').resize((int(w*resc),int(z*resc*aspect)),
                         resample=resample),
                         box=(0,int(h*resc+resc)))

        #add STORM points
        draw=PIL.ImageDraw.Draw(im)
        r=int(stormsize/res)

        if s.nlp>=1:
            self.slices_storm_xy=numpy.empty((s.nlp,2))
            self.slices_storm_xz=numpy.empty((s.nlp,2))
            self.slices_storm_zy=numpy.empty((s.nlp,2))
            cfint,hashes=[],{}
            for i in s.l:
                # STORM's x and y coordinates (normal directions, as opposed to the image's inverted directions
                # where x=height, y=width here. Because of handling numpy arrays.)                
                x=int((float(s.xc(i)-self.offset[0])/self.calib-xmin)*resc+0.5)
                y=int((float(s.yc(i)-self.offset[1])/self.calib-ymin)*resc+0.5)                
                self.slices_storm_xy[i,:]=[x,y]         
                x_xz=int((w+1)*resc+aspect*(midslice-zmipstart)*resc+aspect*(s.zc(i)+zoffset)*resc/zcalib+0.5)                
                self.slices_storm_zy[i,:]=[x_xz,y]                               
                y_zy=int((h+1)*resc+aspect*(midslice-zmipstart)*resc+aspect*(s.zc(i)+zoffset)*resc/zcalib+0.5)                
                self.slices_storm_xz[i,:]=[x,y_zy]
                if not hidestorm:
                    if stormtype=='outline':
                        draw.ellipse((x-r, y-r, x+r, y+r),
                                     outline=colors[s.chnum[i]])
                        draw.ellipse((x_xz, y, x_xz+r, y+r),
                                     outline=colors[s.chnum[i]])
                        draw.ellipse((x-r, y_zy-r, x+r, y_zy+r),
                                     outline=colors[s.chnum[i]])
                    if stormtype=='fill':
                        draw.ellipse((x-r, y-r, x+r, y+r),
                                     fill=colors[s.chnum[i]])
                        draw.ellipse((x_xz-r, y-r, x_xz+r, y+r),
                                     fill=colors[s.chnum[i]])
                        draw.ellipse((x-r, y_zy-r, x+r, y_zy+r),
                                     fill=colors[s.chnum[i]])
                if not measure is None:
                    x=int(float(s.xc(i)-self.offset[0])/self.calib-xmin)
                    y=int(float(s.yc(i)-self.offset[1])/self.calib-ymin)
                    z=int((midslice-zmipstart)+(s.zc(i)+zoffset)/zcalib)
                    hsh=x+y*xmax+z*xmax*ymax
                    if hashes.has_key(hsh):
                        if not (x,y,z) in hashes[hsh]:
                            cfint.append(roi[z,measure,y,x])
                            hashes[hsh].append((x,y,z))
                    else:
                        cfint.append(roi[z,measure,y,x])
                        hashes[hsh]=[(x,y,z)]
            self.cfint_slices=cfint
                        
        #add roi points
        if not hideroi:
            pg=[] 
            for i in xrange(len(self.roi[0])):
                x=int((float(self.roi[0][i]-self.offset[0])/self.calib-xmin)*resc+0.5)
                y=int((float(self.roi[1][i]-self.offset[1])/self.calib-ymin)*resc+0.5)
                im.putpixel((x,y), (255,255,255))
                draw.ellipse((x-3, y-3, x+3, y+3),fill=(255,170,42))
                pg.append((x,y))
            if poly:
                draw.polygon(pg)
        if not hidelasso:
            lasso_numbers=len(self.lasso)
            for n in range(lasso_numbers):
                lasso_poly=[]
                for i in self.lasso[n]:
                    lasso_poly.append((i[0],i[1]))
                draw.polygon(lasso_poly)   
        self.im=im
    
    def selector2D(self,wintitle="Select ROIs", zoom=False):
#    'Select ROIs'):   
        '''Displays the 2D view generated by the showstorm function,
        and allows hand selection of points. Indices of points included
        in each selection are returned.'''
        incl=[]
        #start drawing gui
        fig=pylab.figure()
        pylab.imshow(self.im, origin='upper',aspect='equal')
        ax=fig.add_subplot(111)
        cursor=drawROI(ax,fig)
        fig.canvas.mpl_connect('motion_notify_event', cursor.motion_notify_callback)
        fig.canvas.mpl_connect('button_press_event', cursor.button_press_callback)
        fig.canvas.mpl_connect('key_press_event', cursor.key_press_callback)
        plt.get_current_fig_manager().window.showMaximized()
        fig.canvas.set_window_title(wintitle)
#        pylab.show(block=True)
        global spyder_blocking_bug
        spyder_blocking_bug=True
        while not globals()['spyder_blocking_bug']==False:
            plt.pause(0.5)
        self.selectortrash=cursor.trash
        filtrois = []   
        for points in cursor.lassocoordset:
            thisfilt=[]
            filtrois.append(thisfilt)
            r=numpy.array(points)
            poly=mplpath.Path(r)
            storm_pix=self.showstorm_storm_xy
            for i in self.storm.l:
                if poly.contains_point(storm_pix[i,:]):
                    thisfilt.append(i)                                                      
        for roi in filtrois:
            for point in roi:
                incl.append(point)
        incl=sorted(incl)   #if mulitple rois are selected. 1 point only can be in 1 roi
        plt.close(fig)                
        self.selector2D_incl=incl
        self.lasso=cursor.lassocoordset
        return incl
    
    def selector3D(self, wintitle='Select STORM points to include', ch=None, zoom=False):
        '''Displays the slices view generated by the slices function,
        and allows hand selection of points. Indices of points included
        in all selections are returned.
        If no selection is made, all points are returned.'''
        if ch==None:
            ch=self.storm.chdict
        incl=[]
        #start drawing gui
        fig=pylab.figure()
        pylab.imshow(self.im, origin='upper',aspect='equal')
        ax=fig.add_subplot(111)
        cursor=drawROI(ax,fig)
        fig.canvas.mpl_connect('motion_notify_event', cursor.motion_notify_callback)
        fig.canvas.mpl_connect('button_press_event', cursor.button_press_callback)
        fig.canvas.mpl_connect('key_press_event', cursor.key_press_callback)
        if zoom:
            fig.canvas.mpl_connect('scroll_event',cursor.zoom)
        plt.get_current_fig_manager().window.showMaximized()
        fig.canvas.set_window_title(wintitle)
#        pylab.show(block=True)
        global spyder_blocking_bug
        spyder_blocking_bug=True
        while not globals()['spyder_blocking_bug']==False:
            plt.pause(0.5)
        #find point inside rois or return all if no roi was drawn
        self.selectortrash=cursor.trash
        if len(cursor.lassocoordset)>0:
            filtrois=[]
            for points in cursor.lassocoordset:
                thisfilt=[]
                filtrois.append(thisfilt)
                r=numpy.array(points)
                poly=mplpath.Path(r)
                #detect which view
                if r[0,0]>self.slices_start_zy:
                    storm_pix=self.slices_storm_zy
                elif r[0,1]>self.slices_start_xz:
                    storm_pix=self.slices_storm_xz
                else:
                    storm_pix=self.slices_storm_xy
                for i in self.storm.l:
                    if poly.contains_point(storm_pix[i,:]):
                        thisfilt.append(i)
            #find points which present in all rois
            for i in self.storm.l:
                if not self.storm.channels[i] in ch:
                    continue
                n=0
                for j in filtrois:
                    for k in j:
                        if i==k:
                            n+=1
                if n==len(filtrois):
                    incl.append(i)
            self.selector3D_filtrois=filtrois
        else:
            incl=self.storm.l
        plt.close(fig)        
        self.selector3D_incl=incl
        self.lasso=cursor.lassocoordset
        return incl
    
    def pixselector3D(self, image=None, res=50, wintitle='Select area to include'):
        '''Allows area selection on the 3d image generated with the
        slices function. Returns a coords object with the included points
        Parameter res specifies the spacing of included points.'''
        if image is None:
            self.slices(hidestorm=True)
            image=self.im
        #start drawing gui
        fig=pylab.figure()
        pylab.imshow(image, origin='upper',aspect='equal')
        ax=fig.add_subplot(111)
        cursor=drawROI(ax,fig)
        fig.canvas.mpl_connect('motion_notify_event', cursor.motion_notify_callback)
        fig.canvas.mpl_connect('button_press_event', cursor.button_press_callback)
        fig.canvas.mpl_connect('key_press_event', cursor.key_press_callback)
        plt.get_current_fig_manager().window.showMaximized()
        fig.canvas.set_window_title(wintitle)
        pylab.show(block=True)
        self.selectortrash=cursor.trash
        if len(cursor.lassocoordset)>=3: 
        #move rois to storm coordinates
            zypolys,xzpolys,xypolys=[],[],[]
            xmin,xmax,ymin,ymax,zmin,zmax=[50000,0]*3                     
            for points in cursor.lassocoordset:    
                r=numpy.array(points)
                if r[0,0]>self.slices_start_zy:
                    zyroi=numpy.empty(r.shape)
                    for i, p in enumerate(points):
                        z,y=p
                        zc=self.slices_zcalib*(z-(self.slices_w+1)*self.slices_resc-self.slices_aspect*self.slices_midslice*self.slices_resc)/(self.slices_resc*self.slices_aspect)-self.slices_zoffset
                        yc=(float(y)/self.slices_resc+self.slices_ymin)*self.calib+self.offset[1]
                        zyroi[i,:]=[zc,yc]
                        zmin,zmax=min(zmin,zc),max(zmax,zc)
                        ymin,ymax=min(ymin,yc),max(ymax,yc)
                    zypolys.append(mplpath.Path(zyroi))
                elif r[0,1]>self.slices_start_xz:
                    xzroi=numpy.empty(r.shape)
                    for i, p in enumerate(points):
                        x,z=p
                        xc=(float(x)/self.slices_resc+self.slices_xmin)*self.calib+self.offset[0]
                        zc=self.slices_zcalib*(z-(self.slices_h+1)*self.slices_resc-self.slices_aspect*self.slices_midslice*self.slices_resc)/(self.slices_resc*self.slices_aspect)-self.slices_zoffset
                        xzroi[i,:]=[xc,zc]
                        zmin,zmax=min(zmin,zc),max(zmax,zc)
                        xmin,xmax=min(xmin,xc),max(xmax,xc)
                    xzpolys.append(mplpath.Path(xzroi))
                else:
                    xyroi=numpy.empty(r.shape)
                    for i, p in enumerate(points):
                        x,y=p
                        xc=(float(x)/self.slices_resc+self.slices_xmin)*self.calib+self.offset[0]
                        yc=(float(y)/self.slices_resc+self.slices_ymin)*self.calib+self.offset[1]
                        xyroi[i,:]=[xc,yc]
                        xmin,xmax=min(xmin,xc),max(xmax,xc)
                        ymin,ymax=min(ymin,yc),max(ymax,yc)
                    xypolys.append(mplpath.Path(xyroi))
            xmin,xmax,ymin,ymax,zmin,zmax=map(int,[xmin,xmax,ymin,ymax,zmin,zmax])
            #pick points inside selection
            picked=[]
            for x in xrange(xmin,xmax+res,res):
                for y in xrange(ymin,ymax+res,res):
                    a=False
                    for poly in xypolys:
                        a=poly.contains_point([x,y])
                        if a:
                            break
                    if not a:
                        continue
                    for z in xrange(zmin,zmax+res,res):
                        b=False
                        for poly in xzpolys:
                            b=poly.contains_point([x,z])
                            if b:
                                break
                        if not b:
                            continue
                        c=False
                        for poly in zypolys:
                            c=poly.contains_point([z,y])
                            if c:
                                break
                        if c:
                            picked.append([x,y,z])
                        
            #return data
            if len(picked)>0:
                newstorm=coords()
                newstorm.fill(numpy.asarray(picked).transpose())
                self.pixselector3D_storm=newstorm
                return newstorm
        
    def roifilt(self, storm):
        '''finds storm points within roi as coords class'''
        l=len(self.roi[0])
        r=numpy.empty((l,2))
        for i in xrange(l):
            r[i,:]=[self.roi[0][i],self.roi[1][i]]
        self.poly=mplpath.Path(r)
        incl=[]
        for i in storm.l:
            if self.poly.contains_point(storm.getc(i)[:2]):
                incl.append(i)
        self.filtstorm=storm.getsub(incl)
        return self.filtstorm
        
    def roicfint(self):
        '''measures confocal intensity within ROI'''
        #find bounds
        xmin=int((min(self.roi[0])-self.offset[0])/self.calib)
        xmax=int((max(self.roi[0])-self.offset[0])/self.calib)
        ymin=int((min(self.roi[1])-self.offset[1])/self.calib)
        ymax=int((max(self.roi[1])-self.offset[1])/self.calib)
        w=int(xmax)-int(xmin)+1
        h=int(ymax)-int(ymin)+1
        #make polygon of ROI
        l=len(self.roi[0])
        r=numpy.empty((l,2))
        for i in xrange(l):
            r[i,:]=[self.roi[0][i],self.roi[1][i]]
        self.poly=mplpath.Path(r)
        #visit pixels
        chs=self.cf.ConfocalMetaData['ChannelNum']
        im=numpy.zeros((h,w,chs), dtype=numpy.uint8)
        pixn=0
        for x in xrange(xmin,xmin+w):
            for y in xrange(ymin,ymin+h):
                if self.poly.contains_point([(x+0.5)*self.calib+self.offset[0],(y+0.5)*self.calib+self.offset[1]]):
                    im[y-ymin,x-xmin,:]=self.cf.ConfocalData[:,y,x]
                    pixn+=1
        self.roi_ints=im
        self.roim=PIL.Image.fromarray(im, 'RGB')
        return pixn
        
class drawROI:
    
    def __init__(self, ax, fig):
        self.previous_point = []
        self.start_point = []
        self.end_point = []
        self.line = None    
        self.lassocoords=[]
        self.lassocoordset=[]
        self.prevline=[]
        self.fig =  fig
        self.ax=ax
        self.fig.canvas.draw()
        self.trash=False
        self.orig_xlim=self.ax.get_xlim()
        self.orig_ylim=self.ax.get_ylim()
        
    def key_press_callback(self, event):
        if event.key=='z':
            if len(self.lassocoordset)>0:
                self.lassocoordset.pop(-1)
                line=self.prevline.pop(-1)
                line.set_color('r')
                event.inaxes.add_line(line)              
                self.fig.canvas.draw()
        elif event.key=='q':
            global spyder_blocking_bug
            globals()['spyder_blocking_bug']=False
            cfg.spyder_blocking_bug = False
#            plt.close(event.canvas.figure)
        elif event.key=='x':
#            self.trash=True
            global spyder_blocking_bug
            globals()['spyder_blocking_bug']=False
            cfg.spyder_blocking_bug = False
#            plt.close(event.canvas.figure)
        elif event.key=='r':
            self.ax.set_xlim(self.orig_xlim)
            self.ax.set_ylim(self.orig_ylim)
        
    def motion_notify_callback(self, event):
        if event.inaxes:
            ax = event.inaxes
            x, y = event.xdata, event.ydata
#            if event.key=="0": #enable line if drawing should only occur when pressing 0
            if event.button == None and self.line != None: # Move line around 
                self.line.set_data([self.previous_point[0], x],
                                   [self.previous_point[1], y]) 
                self.fig.canvas.draw()
                
            elif event.button == 1: # Free Hand Drawing
                line = pylab.Line2D([self.previous_point[0], x],
                              [self.previous_point[1], y])                  
                self.lassocoords.append([x,y])
                ax.add_line(line)
                self.previous_point = [x, y]
                self.fig.canvas.draw()

        
    def button_press_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            ax = event.inaxes 
#            if event.key=="0":    #enable line if drawing should only occur when pressing 0            
            if event.button == 1:  # If you press the left button
                if self.line == None: # if there is no line, create a line
                    self.line = pylab.Line2D([x,  x],
                                       [y, y],
                                       marker = 'o')
                    self.lassocoords.append([x,y])
                    self.start_point = [x,y]
                    self.previous_point =  self.start_point 
                    ax.add_line(self.line)
                    self.fig.canvas.draw()
                # add a segment
                else: # if there is a line, create a segment
                    self.line = pylab.Line2D([self.previous_point[0], x], 
                                             [self.previous_point[1], y],
                                             marker = 'o')
                    self.lassocoords.append([x,y])
                    self.previous_point = [x,y]
                    event.inaxes.add_line(self.line)
                    self.fig.canvas.draw()
            
            elif event.button == 3 and self.line != None: # close the loop
                        self.line.set_data([self.previous_point[0], self.start_point[0]],
                                           [self.previous_point[1], self.start_point[1]])                       
                        
                        self.fig.canvas.draw()
                        self.prevline.append(pylab.Line2D([self.previous_point[0],self.start_point[0]],[self.previous_point[1],self.start_point[1]]))
                        self.line = None
                        ## self.lassocoords.append([x,y]) 
                        ## no need to append this! Last point for closing (previous_point) already 
                        ## in lassocoord; otherwise, a point where mouse right click was clicked 
                        ## will be added too, N points in lasso will be 1 bigger than vertices
                        ## you wanted to add
                        self.lassocoordset.append(copy.deepcopy(self.lassocoords))
                        self.lassocoords=[]
                        self.fig.canvas.setFocus() 

#    def zoom(self, event): # https://stackoverflow.com/questions/11551049/matplotlib-plot-zooming-with-scroll-wheel
#            if event.inaxes:
#                base_scale=1.2
#                cur_xlim = self.ax.get_xlim()
#                cur_ylim = self.ax.get_ylim()
#                cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
#                cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
#                xdata = event.xdata 
#                ydata = event.ydata 
#                if event.button == 'up':
#                    scale_factor = 1/base_scale
#                elif event.button == 'down':
#                    scale_factor = base_scale
#                else:
#                    scale_factor = 1
#                    print event.button
#                self.ax.set_xlim([xdata - cur_xrange*scale_factor,
#                             xdata + cur_xrange*scale_factor])
#                self.ax.set_ylim([ydata - cur_yrange*scale_factor,
#                             ydata + cur_yrange*scale_factor])
#                plt.draw()
        
def dist(c1,c2):
    '''distance of 2 3d points provided as coordinate triplets'''
    x,y,z=map(operator.sub,c1,c2)
    return math.sqrt(x*x+y*y+z*z)