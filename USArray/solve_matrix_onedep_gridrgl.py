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
    print(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(dlon))
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

    # change path
    os.chdir('{}'.format(dirname))
    
    # "10s" means a single 10-byte string, "10c" means 10 characters
    bihead_format = "i" + "8s" + "f"*6

    target_depth = 450 # km
    with open("{}".format(filename),"rb") as fp:
        while True:
            bihead = fp.read(36)
            if bihead:
                head = struct.unpack(bihead_format,bihead)

                if ( np.abs(np.array(head[2])-target_depth) <0.1 and head[6]<2.0 and np.abs(head[5])<1):
                    grid_j.append(head[0])
                    a = head[1].decode()
                    evtcd.append(a)
                    grid_lon.append(head[3])
                    grid_lat.append(head[4])
                    results.append(head[5])
                    gcarc.append(head[6])
#                    print(head[0],head[2],head[3],head[4],head[5])
            else:
                break
    grid_j = np.array(grid_j)
    evtcd = np.array(evtcd)
    grid_lon = np.array(grid_lon)
    grid_lat = np.array(grid_lat)
    results = np.array(results)

    print(len(evt))
    evt = evt[seismnum>minseismnum]
    print(len(evt))
    number_event = len(evt)

    # Regularization of the grid
    lamda0 = 10.0
    delta0 = 1.0
    col_reg = []
    row_reg = []
    value_reg = []
    m = 0
    for i in np.arange(number_grid):
        for j in np.arange(i,number_grid):
            deltaijkl = degdist(gridpoints[i,1],gridpoints[i,0],gridpoints[j,1],gridpoints[j,0])
            if (deltaijkl < 5 and deltaijkl > 0.000001):
                row_reg.append(m)
                col_reg.append(i)
                value_reg.append(lamda0*np.exp(-deltaijkl/delta0))
                row_reg.append(m)
                col_reg.append(j)
                value_reg.append(-lamda0*np.exp(-deltaijkl/delta0))
                m = m + 1
    # Regularization of events
    evtall = pd.read_csv('{}/events_cat_all.txt'.format(savepath),sep=' ',header=None,skipinitialspace=True)
    evtall_code = evtall[0]
    evtall_lat = evtall[2]
    evtall_lon = evtall[3]
    for i in np.arange(number_event):
        for j in np.arange(i,number_event):
            idx1 = list(evtall_code).index(evt[i])
            idx2 = list(evtall_code).index(evt[j])
            deltaijkl = degdist(evtall_lat[idx1],evtall_lon[idx1],evtall_lat[idx2],evtall_lon[idx2])
            if (deltaijkl < 5 and deltaijkl > 0.000001):
                row_reg.append(m)
                col_reg.append(i+number_event)
                value_reg.append(lamda0*np.exp(-deltaijkl/delta0))
                row_reg.append(m)
                col_reg.append(j+number_event)
                value_reg.append(-lamda0*np.exp(-deltaijkl/delta0))
                m = m + 1
    value_reg = np.array(value_reg)
    A_reg = coo_matrix((value_reg,(row_reg,col_reg)),shape=(m,number_grid+number_event))
    rhs_reg = np.zeros(m)

    # Measurements matrix
    X = []
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
    print(rhs0.shape)
    A0 = coo_matrix((np.ones(len(row)),(row,col)),shape=(len(rhs0),number_grid+number_event))

    
    # Least square problem
    A = vstack((A0,A_reg))
    rhs = np.concatenate((rhs0,rhs_reg),axis=0)
    x,istop,itn,r1norm,r2norm,anorm,acond,arnorm,xnorm,var = lsqr(A,rhs,atol=1e-12,btol=1e-12)
    x[np.where(count<1)] = 0
    X.append(x[0:number_grid])
    X = np.array(X)
    print(x,np.amax(X))
    print(np.linalg.norm(rhs[0:number_grid+number_event]))
    print(np.linalg.norm((A*x-rhs)[0:number_grid+number_event]))
#    np.savetxt('{}/matrix_rhs.txt'.format(dirname),rhs,fmt="%8.3f")
#    np.savetxt('{}/matrix_A.txt'.format(dirname),coo_matrix.todense(A),fmt="%8.3f")

    # plot results
    fig,ax = plt.subplots(1,1,figsize=[8,4])
    cm = ax.imshow(X.reshape(ncol,nrow).T,cmap='bwr',vmin=-0.05,vmax=0.05,aspect='auto',origin='lower')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_yticks([0,5,10,15,20])
    ax.set_yticklabels(['25','30','35','40','45','50'])
    ax.set_xticks([2,12,22,32,42,52])
    ax.set_xticklabels(['-125','-115','-105','-95','-85','-75'])
    ax.set_title('Reflection depth = {} km'.format(target_depth))
    plt.colorbar(cm,ax=ax)
    plt.savefig('{}/solve_matrix_onedep_gridrgl.pdf'.format(savepath))
    plt.close()

    # plot count map
    fig,ax = plt.subplots(1,1,figsize=[8,4])
    cm = ax.imshow(count.reshape(ncol,nrow).T,cmap='rainbow',aspect='auto',origin='low',norm=colors.LogNorm())
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_yticks([0,5,10,15,20])
    ax.set_yticklabels(['25','30','35','40','45','50'])
    ax.set_xticks([2,12,22,32,42,52])
    ax.set_xticklabels(['-125','-115','-105','-95','-85','-75'])
    plt.colorbar(cm,ax=ax)
    plt.savefig('{}/solve_matrix_onedep_count.pdf'.format(savepath))
    plt.close()

    # save event results
    evtout = []
    for i in np.arange(number_event):
        idx = list(evtall_code).index(evt[i])
        evtout.append([evtall_code[idx],evtall_lat[idx], evtall_lon[idx], x[i+number_grid]])
    evtout = np.array(evtout)
    np.savetxt('{}/evtout.txt'.format(savepath),evtout,fmt="%s %s %s %s")

    # save gridpoint results
    gridpoint_terms = []
    for i in np.arange(number_grid):
        gridpoint_terms.append([gridpoints[i,0],gridpoints[i,1],x[i]])
    # save as latitude, longitude, inversion results
    np.savetxt('{}/gridpointout.txt'.format(savepath),gridpoint_terms,fmt="%8.3f %8.3f %8.3f")


def main():
    filename = 'rp_Asia_interval1_area_deg2_win10_150_800.BHT.vel.snrS.xh'
    xhfilename = 'all_Asia.BHT.vel.snrS.xh'
    dirname = '/home/meichen/work1/SH_TOMO/Combined_xh'
    curpath = '/home/meichen/Research/SH_TOMO/USArray'
    minseismnum = 10

    evt,seismnum = readxh(dirname=dirname,filename=xhfilename)
    readbinary(filename=filename,dirname=dirname,savepath=curpath,evt=evt,seismnum=seismnum,minseismnum=minseismnum)

main()
