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
    mincount = kwargs.get('mincount')

    gridpoints = np.genfromtxt('{}/grid.txt'.format(dirname))
    number_grid = len(gridpoints)
    nrow = len(np.unique(gridpoints[:,1]))
    ncol = len(np.unique(gridpoints[:,0]))
    evt = evt[seismnum>minseismnum]
    number_event = len(evt)

    # Regularization matrix
    lamda0 = 10.0
    delta0 = 1.0
    col_reg = []
    row_reg = []
    value_reg = []
    m = 0
    for i in np.arange(number_grid):
        for j in np.arange(i,number_grid):
            deltaijkl = degdist(gridpoints[i,1],gridpoints[i,0],gridpoints[j,1],gridpoints[j,0])
            if (deltaijkl < 4 and deltaijkl > 0.000001):
                row_reg.append(m)
                col_reg.append(i)
                value_reg.append(lamda0*np.exp(-deltaijkl/delta0))
                row_reg.append(m)
                col_reg.append(j)
                value_reg.append(-lamda0*np.exp(-deltaijkl/delta0))
                m = m + 1
    value_reg = np.array(value_reg)
    rhs_reg = np.zeros(m)

    os.chdir('{}'.format(dirname))    

    # "10s" means a single 10-byte string, "10c" means 10 characters
    bihead_format = "i" + "8s" + "f"*7

    X = np.zeros((199,number_grid))
    #dep = np.hstack((np.arange(355,455,5),np.arange(605,755,5)))
    dep = np.arange(10,1005,5)
    for target_depth in dep:
        idep = (int) (target_depth/5-2)
        grid_j = []
        evtcd = []
        grid_lon = []
        grid_lat = []
        results = []
        gcarc = []
        az = []
        evtstngcarc = []
        with open("{}".format(filename),"rb") as fp:
            while True:
                bihead = fp.read(40)
                if bihead:
                    head = struct.unpack(bihead_format,bihead)
   
                    if (np.abs(np.array(head[2])-target_depth)<0.1) and head[6] < 2.0 and np.abs(head[5]) < 1:
                        grid_j.append(head[0])
                        a = head[1].decode()
                        evtcd.append(a)
                        grid_lon.append(head[3])
                        grid_lat.append(head[4])
                        results.append(head[5])
                else:
                    break
        grid_j = np.array(grid_j)
        evtcd = np.array(evtcd)
        grid_lon = np.array(grid_lon)
        grid_lat = np.array(grid_lat)
        results = np.array(results)

        lon = np.sort(np.unique(grid_lon))
        lat = np.sort(np.unique(grid_lat))    

        # Regularization matrix
        A_reg = coo_matrix((value_reg,(row_reg,col_reg)),shape=(m,number_grid+number_event))

        # Measurement matrix
        row = []
        col = []
        rhs0 = []
        count = np.zeros(number_grid)
        row_line = 0
        for i in np.arange(len(grid_j)):
            for j in np.arange(number_event):
                if (evtcd[i] == evt[j]):
                    row.append(row_line)
                    col.append(grid_j[i])
                    row.append(row_line)
                    col.append(j+number_grid)
                    rhs0.append(results[i])
                    count[grid_j[i]] = count[grid_j[i]] + 1
                    row_line = row_line + 1
        rhs0 = np.array(rhs0)
        A0 = coo_matrix((np.ones(len(row)),(row,col)),shape=(len(rhs0),number_grid+number_event))

        # Build matrix
        A = vstack((A0,A_reg))
        rhs = np.concatenate((rhs0,rhs_reg),axis=0)
        x = np.zeros(A.shape[1])

        # extract grid points with enough numbers of counts
        extract_index = np.hstack((np.where(count>mincount)[0],np.arange(number_grid,number_grid+number_event)))
        A = A.tocsr()[:,extract_index].tocoo()

        # Solve least square problem
        x[extract_index] = lsqr(A,rhs,atol=1e-12,btol=1e-12)[0]
        X[idep,:] = x[0:number_grid]
        print(target_depth,number_grid,x)
    X = np.array(X)
    np.savetxt('results.xyz',np.hstack((np.arange(10,1005,5).reshape(199,1),X)),fmt="%12.8f")
    X = X.reshape(199,ncol,nrow)

    const_lat_list = np.array([25.0,27.0,29.0,31.0,33.0,35.0,37.0,39.0,41.0,43.0,45.0,47.0,49.0])
    const_lon_list = np.array([-125.0,-120.0,-115.0,-110.0,-105.0,-100.0,-95.0,-90.0,-85.0,-80.0,-75.0])
    fig,ax = plt.subplots(max(len(const_lat_list),len(const_lon_list)),2,figsize=[16,5*max(len(const_lon_list),len(const_lat_list))])
    for i in np.arange(len(const_lat_list)):
        const_lat = const_lat_list[i]
        k = list(np.sort(np.unique(gridpoints[:,1]))).index(const_lat)

        cm = ax[i,0].imshow(X[:,:,k],cmap='bwr',vmin=-0.1,vmax=0.1,aspect='auto')
        ax[i,0].set_xlabel('Longitude')
        ax[i,0].set_ylabel('Depth')
        ax[i,0].set_title('lat={}'.format(const_lat))
        #ax.set_yticks([0,5,10,15,20])
        #ax.set_yticklabels(['25','30','35','40','45','50'])
        #ax.set_xticks([2,12,22,32,42,52])
        #ax.set_xticklabels(['-125','-115','-105','-95','-85','-75'])
    for i in np.arange(len(const_lon_list)):
        const_lon = const_lon_list[i]
        k = list(np.sort(np.unique(gridpoints[:,0]))).index(const_lon)

        cm = ax[i,1].imshow(X[:,k,:],cmap='bwr',vmin=-0.1,vmax=0.1,aspect='auto')
        ax[i,1].set_xlabel('Latitude')
        ax[i,1].set_ylabel('Depth')
        ax[i,1].set_title('lon={}'.format(const_lon))
    plt.colorbar(cm,ax=ax[0,0])
    plt.savefig('{}/{}.pdf'.format(savepath,filename))
    plt.close()


def main():
    filename = 'rp_s40rts_10s_interval1_area_deg3_win5_10_1000_dep5_deg_50_110.MXT.vel.snrS.xh'
    xhfilename = 's40rts_10s.snrS.MXT.vel.xh'
    dirname = '/home/meichen/work1/SH_TOMO/synthetics/XH_files'
    curpath = '/home/meichen/Research/SH_TOMO/synthetics'
    minseismnum = 10
    mincount=5

    evt,seismnum = readxh(dirname=dirname,filename=xhfilename)
    readbinary(filename=filename,dirname=dirname,savepath=curpath,evt=evt,seismnum=seismnum,minseismnum=minseismnum,mincount=mincount)

main()
