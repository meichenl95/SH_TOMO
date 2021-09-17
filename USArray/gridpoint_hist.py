#!/home/meichen/anaconda3/bin/python

def root_stack(seism,n):

    import numpy as np

    seism_root = np.zeros(len(seism))
    for i in np.arange(len(seism)):
        sign = seism[i]/np.abs(seism[i])
        seism_root[i] = sign*np.sqrt(np.abs(seism[i]))
    return seism_root

def degdist(lat1,lon1,lat2,lon2):

## latitude and longitude are in degrees

    import math
    import numpy as np

    lat1 = np.deg2rad(lat1)
    lon1 = np.deg2rad(lon1)
    lat2 = np.deg2rad(lat2)
    lon2 = np.deg2rad(lon2)
    dlon = lon2 - lon1

    deg = np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(dlon))
    deg = np.abs(np.rad2deg(deg))

    return deg

def readbinary(**kwargs):

# parameters
# filename		the name of the binary file to read in
# dirname		the directory of the binary file
# savepath		the directory to save output
# ndata			the number of data in each seismogra
# N1			the starting point of the window interested in
# N2			the ending point of the window interested in
# stack_num		the minimum number of stacked seismograms of each grid point

    import struct
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.colors as colors
    import obspy.signal.filter
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import lsqr
    from scipy.sparse import vstack
    import pandas as pd

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    savepath = kwargs.get('savepath')
    ndata = kwargs.get('ndata')
    target_lat = kwargs.get('target_lat')
    target_lon = kwargs.get('target_lon')
    target_depth = kwargs.get('target_depth') # km

    evtcd = []
    results = []

    # change path
    os.chdir('{}'.format(dirname))
    
    # "10s" means a single 10-byte string, "10c" means 10 characters
    bihead_format = "i" + "8s" + "f"*7

    with open("{}".format(filename),"rb") as fp:
        while True:
            bihead = fp.read(40)
            if bihead:
                head = struct.unpack(bihead_format,bihead)

                if ( np.abs(np.array(head[2])-target_depth) <0.1 and head[6]<2.0 and np.abs(head[4]-target_lat)<0.01 and np.abs(head[3]-target_lon)<0.01 and head[8]>80 and head[8]<100):
                    results.append(head[5])
            else:
                break
    results = np.array(results)

    # plot results
    fig,ax = plt.subplots(1,1,figsize=[8,4])
    ax.hist(results,bins=20)
    ax.set_xlabel('Measurements')
    ax.set_ylabel('Counts')
    ax.set_title('lat={}, lon={}, reflect_dep={} km\n mean={:5.3f}, median={:5.3f}, std={:5.3f}'.format(target_lat,target_lon,target_depth,np.mean(results),np.median(results),np.std(results)))
    plt.savefig('{}/gridpoint_hist.pdf'.format(savepath))
    plt.close()

def main():
    filename = 'rp_Asia_interval1_area_deg3_win10_405_450_dep5_deg_50_110.BHT.vel.snrS.xh'
    dirname = '/home/meichen/work1/SH_TOMO/Combined_xh'
    curpath = '/home/meichen/Research/SH_TOMO/USArray'

    readbinary(filename=filename,dirname=dirname,savepath=curpath,target_lat=43,target_lon=-117,target_depth=450)

main()
