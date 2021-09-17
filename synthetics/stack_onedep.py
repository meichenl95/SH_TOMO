#!/home/meichen/anaconda3/bin/python

def readxh(**kwargs):                                                       
 
##---------------------------------------##
 
# This function read the xh file header.
 
# Created by Meichen Liu on July 6th, 2019.

##---------------------------------------##
##parameters
 
 
    import struct
    from obspy.core import UTCDateTime
    import matplotlib.pyplot as plt
    import numpy as np
    import os
 
    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname') 

    os.chdir('{}'.format(dirname))
 
    evtcd = []
    xhhead_format = "fii" + "f"*11 + "iiiiif"*2 + "ifffiii" + "ff"*60 + "fff"+"f"*20*2+"i"*20 + "14s8s8s8s8s72s8s34s"
    with open("{}".format(filename),"rb") as fp:
        while True:
            xhhead = fp.read(1024)
            if xhhead:
                head = struct.unpack(xhhead_format,xhhead)
 
                version = head[0]
                nhdr = head[1] 
                i12345678 = head[2]
                elat = head[3]
                elon = head[4]
                edep = head[5]
                Mb = head[6]
                Ms = head[7]
                Mw = head[8]
                slat = head[9]
                slon = head[10]
                elev = head[11]
                azim = head[12]
                incl = head[13]
                ot = head[14:20]
                tstart = head[20:26]
                                     
                ndata = head[26]
                delta = head[27]
                tshift = head[28]
                maxamp = head[29]
                qual = head[30]
                chid = head[31]
                locc = head[32]
                                 
                poles = head[33:93]
                zeros = head[93:153]
                                     
                DS = head[153]
                A0 = head[154]
                f12345678 = head[155]
                tpck = head[156:176]
                flt = head[176:196]
                intg = head[196:216]
                                 
                cmtcd = head[216]
#                if flt[0]>5 and flt[1]>5 and flt[2]>5 and flt[3]>2.5 and edep<35:
                evtcd.append(''.join(list(head[217].decode())[0:7]))
                netw = head[218]
                stnm = head[219]
                chan = head[220]
                rcomment = head[221] 
                wavf = head[222]
                padding = head[223]
                                     
                xhdata = fp.read(ndata*4)
                xhdata_format = "f" * ndata
                amp = struct.unpack(xhdata_format,xhdata)
            else:
                break
    evt = np.unique(evtcd)
    seismnum = []
    for eventcode in evt:
        recordlist = [i for i,item in enumerate(evtcd) if eventcode == item]
        seismnum.append(len(recordlist))

    evt = np.array(evt)
    seismnum = np.array(seismnum)
    return np.array(evt),np.array(seismnum)

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
    evt = kwargs.get('evt')
    seismnum = kwargs.get('seismnum')
    minseismnum = kwargs.get('minseismnum')
    target_depth = kwargs.get('target_depth')
    average_or_root = kwargs.get('average_or_root')

    gridpoints = np.genfromtxt('{}/grid.txt'.format(dirname))
    number_grid = len(gridpoints)
    nrow = len(np.unique(gridpoints[:,1])) # along latitude
    ncol = len(np.unique(gridpoints[:,0])) # along longitude
    grid_j = []
    evtcd = []
    grid_lon = []
    grid_lat = []
    gcarc = []
    results = []
    az = []
    evtstngcarc = []

    # change path
    os.chdir('{}'.format(dirname))
    
    # "10s" means a single 10-byte string, "10c" means 10 characters
    bihead_format = "i" + "8s" + "f"*12

    evtall = pd.read_csv('~/Research/SH_TOMO/synthetics/events_cat.txt'.format(savepath),sep=' ',header=None,skipinitialspace=True)
    evtall_code = evtall[0]
    evtall_lat = evtall[2]
    evtall_lon = evtall[3]
    evtall_dep = evtall[4]

    with open("{}".format(filename),"rb") as fp:
        while True:
            bihead = fp.read(60)
            if bihead:
                head = struct.unpack(bihead_format,bihead)
                a = head[1].decode()[0:7]
                idx = list(evtall_code).index(a)
                event_dep = evtall_dep[idx]

                if ( np.abs(np.array(head[2])-target_depth) <0.1):
                    grid_j.append(head[0])
                    evtcd.append(a)
                    grid_lon.append(head[3])
                    grid_lat.append(head[4])
                    if (average_or_root == 'average'):
                        results.append(head[5])
                    elif (average_or_root == 'max'):
                        results.append(head[13])
                    gcarc.append(head[6])
                    az.append(head[7])
                    evtstngcarc.append(head[8])
#                    print(head[0],head[2],head[3],head[4],head[5])
            else:
                break
    grid_j = np.array(grid_j)
    evtcd = np.array(evtcd)
    grid_lon = np.array(grid_lon)
    grid_lat = np.array(grid_lat)
    results = np.array(results)
    evtstngcarc = np.array(evtstngcarc)

    unique = np.unique(np.vstack((results,evtstngcarc,evtcd)),axis=1).T
    evt = evt[seismnum>minseismnum]
    number_event = len(evt)

    # Measurements matrix
    event_term = np.zeros(number_event)
    num_event = np.zeros(number_event)
    for i in np.arange(unique.shape[0]):
        for j in np.arange(number_event):
            if (unique[i,2] == evt[j]):
                event_term[j] = event_term[j] + float(unique[i,0])
                num_event[j] = num_event[j] + 1
    idx = np.where(num_event>0)[0]
    event_term[idx] = event_term[idx] / num_event[idx] /2.0
    

    grid_term = np.zeros(number_grid)
    num_grid = np.zeros(number_grid)
    for i in np.arange(len(grid_j)):
        for j in np.arange(number_event):
            if (evtcd[i] == evt[j]):
                grid_term[int(grid_j[i])] = grid_term[int(grid_j[i])] + results[i] - event_term[j]
                num_grid[int(grid_j[i])] = num_grid[int(grid_j[i])] + 1
    idx = np.where(num_grid>0)[0]
    grid_term[idx] = grid_term[idx] / num_grid[idx]

    # smooth grid terms
    smooth_grid = np.zeros(number_grid)
    smooth_num = np.zeros(number_grid)
    for i in np.arange(number_grid):
        for j in np.arange(number_grid):
            if abs(gridpoints[i,0]-gridpoints[j,0])<=2.5 and abs(gridpoints[i,1]-gridpoints[j,1])<=2.5:
                smooth_grid[i] = smooth_grid[i] + grid_term[j]
                smooth_num[i] = smooth_num[i] + 1
    idx = np.where(smooth_num>0)[0]
    smooth_grid[idx] = smooth_grid[idx] / smooth_num[idx]

    # save event results
    evtout = []
    for i in np.arange(number_event):
        try:
            idx = list(evtall_code).index(evt[i])
            evtout.append([evtall_code[idx],evtall_lat[idx], evtall_lon[idx], event_term[i]])
        except:
            continue
    evtout = np.array(evtout)
    np.savetxt('{}/evtout_stack.txt'.format(savepath),evtout,fmt="%s %s %s %s")

    # save gridpoint results
    gridpoint_terms = []
    for i in np.arange(number_grid):
        gridpoint_terms.append([gridpoints[i,0],gridpoints[i,1],smooth_grid[i],num_grid[i]])
    # save as latitude, longitude, inversion results, count number
    np.savetxt('{}/gridpointout_{:04d}_stack_ttcorr.txt'.format(savepath,int(target_depth)),gridpoint_terms,fmt="%8.3f %8.3f %8.5f %d")


def main():

    import sys
    import numpy as np

    xhfilename = 's40rts_nocrust_10s.snrS.BXT.dis.xh'
    dirname = '/home/meichen/work1/SH_TOMO/synthetics/XH_files'
    curpath = '/home/meichen/Research/SH_TOMO/synthetics'
    minseismnum = 0
    average_or_root = 'average'

    evt,seismnum = readxh(dirname='{}/S40RTS_NOCRUST_10s'.format(dirname),filename=xhfilename)

    for target_depth in np.arange(10,1001,5):
        print(target_depth)
        filename = 'dis15.ttcorr.{:04d}'.format(int(target_depth))
        readbinary(filename=filename,dirname='{}/rp_S40RTS_NOCRUST_10s_ttcorr'.format(dirname),savepath=curpath,minseismnum=minseismnum,target_depth=target_depth,average_or_root=average_or_root,evt=evt,seismnum=seismnum)

main()
