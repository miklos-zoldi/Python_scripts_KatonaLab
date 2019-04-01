import sl_mod_miki as sl
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab
from matplotlib.path import Path
import cfg
### %matplotlib auto at the console

def get_middle_coords(filename, center_zoom_factor=0.05 ):
    '''read storm coordinate file and crop data by center_zoom_factor: 
    Returns:
    1: whole cropped dataset as pandas dataframe
    2: cropped [[Xwc, Ywc, Zc]] values as numpy arrays
    '''
    
    #read data as pandas dataframes and get xyz values as numpy arrays 
    df = pd.read_csv(filename, delimiter='\t')
    df = df[df['Channel Name'].astype(str) != 'Z Rejected'] 
    xyz = df[['Xwc','Ywc', 'Zc']].values
    #crop specified area of image
    if center_zoom_factor:
        box = np.array([np.min(xyz, axis=0),
                         np.max(xyz, axis=0)])
        box_center = 0.5 * (box[0, :] + box[1, :])
        box_edges = box[1, :] - box[0, :]
        factor = center_zoom_factor * 0.5
        crop_box = np.array([box_center - factor * box_edges,
                             box_center + factor * box_edges])
        greater_mask = np.all(crop_box[0, :2] < xyz[:, :2], axis=1)
        less_mask = np.all(xyz[:, :2] < crop_box[1, :2], axis=1)
        mask = greater_mask * less_mask
        xyz = xyz[mask, :]
        df = df[mask]
    return df, xyz

basepath = 'L:/Miki/STORM_projects/190301_z-stack_mitochondrium//' 
path = basepath + 'mito_ROIs//'
npath = path #basepath + 'refine_manual_try//'

#def handle_close(event):
#    print " in handle_close_event"
#    global spyder_blocking_bug
#    globals()['spyder_blocking_bug']=False
#    plt.close(event.canvas.figure)

for filename in os.listdir(path)[0:3]:
    if not "RoiCoords" in filename:
        continue
    print filename
    dt, xyz =  get_middle_coords(path+filename, center_zoom_factor=1)
##    b=sl.bouton(path,f,cf=None)
##    b.storm=sl.coords(path,f, rich=True)
###    b.slices()
##    inclpoints=b.selector3D()
##    if not b.selectortrash:
##        b.storm = b.storm.getsub(inclpoints) 
##        b.slices(stormsize=25, stormtype='fill', hideroi=True, hidelasso=False)
##        b.storm.save(npath,f.replace('.tif', '3D.txt'))
###        b.im.save(npath+f.replace('.tif','3D.tif'))
#    
    fig, (ax, ax2) = plt.subplots(nrows=1, ncols=2) #, figsize=(10,10),)
    plt.get_current_fig_manager().window.showMaximized() #max window size
#    plt.ion()
    scatter = ax.scatter(xyz[:,0], xyz[:,1], s=10, c='purple') #xy
    ax.set_aspect(1.0)
    scatter2 = ax2.scatter(xyz[:,0], xyz[:,2], s=10, c='purple') #xz
    ax2.set_aspect(1.0)
#    scatter3 = ax3.scatter(xyz[:,1], xyz[:,2], s=10, c='purple') #yz
    xmin, xmax = xyz[:,0].min(), xyz[:,0].max()
    ymin, ymax = xyz[:,1].min(), xyz[:,1].max()
#    plt.xlim((xmin-500, xmin+500))
#    plt.ylim((ymin-500, ymax+500))

    #instantiate sl drawROI class and get callbacks
    cursor=sl.drawROI(ax,fig)
    fig.canvas.mpl_connect('motion_notify_event', cursor.motion_notify_callback)
    fig.canvas.mpl_connect('button_press_event', cursor.button_press_callback)
    fig.canvas.mpl_connect('key_press_event', cursor.key_press_callback)
#    fig.canvas.mpl_connect('close_event', handle_close)
    

#    while True: #script dies if closing figure!!!
#        if plt.waitforbuttonpress(): #whatetever bouton you press!!!  
#            break
#    global spyder_blocking_bug #pauses figure until spyder_blocking_bug set to False; by default, blocking does not work in Spyder
#    spyder_blocking_bug = True
    while not cfg.spyder_blocking_bug == False:
#        print cfg.spyder_blocking_bug
        plt.pause(0.5)
    #output of freehand selection
    print "got out from loop"
    points_for_polygon = cursor.lassocoordset
    print points_for_polygon
#    print "have points"
#    if len(points_for_polygon) == 0: #No ROI was drawn, continue loop; plot closes without saving
#        print "no ROI was drawn, not saving file, continuing cycle"
#        plt.close()
#        continue
#    
#    #iterate over all drawn polygons
#    x = xyz[:, 0] #2d array for efficient concatenation
#    y = xyz[:, 1]
#    z = xyz[:, 2] 
#    mask = np.ones(len(x), dtype=bool)
#    
#    for poly,view,ax, Nr in zip(points_for_polygon, [y,z], [ax,ax2], [1,2]): #when looping on polys, even if 1 was drawn, no need for squeeze
#        curr_poly = np.array(poly)
##        curr_poly = np.squeeze(curr_poly) #now works with 1 lasso
#        polygon = Path(curr_poly)
#        contain = polygon.contains_points(
#                np.concatenate((x[:, np.newaxis],view[:, np.newaxis]), axis=1))
#        #redraw filtered moleculeset
#        filtered_xyz = xyz[contain]
#        nonFiltered_xyz = xyz[np.logical_not(contain)]
#        ax.scatter(nonFiltered_xyz[:,0], nonFiltered_xyz[:,Nr], s=10, c='white', alpha=0.7)
#        fig.canvas.draw()
#        #update mask
#        mask = mask*contain
#        
#    if not mask.any(): #no points inside any of the ROIs; plot closes without saving
#        print "no point inside ROI, not saving file, continuing cycle"
#        plt.close()
#        continue
#    #filtering for points stayed in all views
#
#    while True: #script dies if closing figure!!!
#        if plt.waitforbuttonpress():
#            plt.close()
#            break
#    filtered_dt = dt[mask]
##    filtered_dt.to_csv(npath+filename.replace('.txt','_Refined.txt'), sep = '\t')

    

    
    
    