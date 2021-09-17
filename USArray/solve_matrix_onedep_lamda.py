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

    target_depth = 400 # km
    with open("{}".format(filename),"rb") as fp:
        while True:
            bihead = fp.read(36)
            if bihead:
                head = struct.unpack(bihead_format,bihead)

                if ( np.abs(np.array(head[2])-target_depth) <0.1 and head[6]<0.5 and np.abs(head[5])<1):
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

    evt = np.unique(evtcd)
    number_event = len(np.unique(evtcd))


    # Measurements matrix
    row = []
    col = []
    rhs0 = []
    count = np.zeros(number_grid)
    for i in np.arange(len(grid_j)):
        for j in np.arange(number_event):
            if (evtcd[i] == evt[j]):
                row.append(i)
                col.append(grid_j[i])
                row.append(i)
                col.append(j+number_grid)
                rhs0.append(results[i])
        count[grid_j[i]] = count[grid_j[i]] + 1
    rhs0 = np.array(rhs0)
#    print(len(grid_j),target_depth)
    A0 = coo_matrix((np.ones(len(row)),(row,col)),shape=(len(grid_j),number_grid+number_event))

    lamda0_list = np.array([1,2,4,8,10,14,20,30,40,50])
    #lamda0_list = np.array([0.1,0.2,0.3,0.4,0.5,0.8,1,2,4,8,12,20,30,40])
    fig,ax = plt.subplots(1+len(lamda0_list),1,figsize=[4,2.5*len(lamda0_list)+2.5])
    misfitreduce = []
    # vary lamda
    for index,lamda0 in enumerate(lamda0_list):
        # Regularization matrix
        delta0 = 1.0
        col_reg = []
        row_reg = []
        value_reg = []
        X = []
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
        value_reg = np.array(value_reg)
        A_reg = coo_matrix((value_reg,(row_reg,col_reg)),shape=(m,number_grid+number_event))
        rhs_reg = np.zeros(m)
    
        # Least square problem
        A = vstack((A0,A_reg))
        rhs = np.concatenate((rhs0,rhs_reg),axis=0)
        x,istop,itn,r1norm,r2norm,anorm,acond,arnorm,xnorm,var = lsqr(A,rhs,atol=1e-12,btol=1e-12)
        X.append(x[0:number_grid])
        X = np.array(X)
        normrhs=np.linalg.norm(rhs)
        misfit=np.linalg.norm((A*x-rhs))
        misfitreduce.append((normrhs-misfit)/normrhs)
        cm = ax[index].imshow(X.reshape(ncol,nrow).T,cmap='bwr',vmin=-0.1,vmax=0.1,aspect='auto',origin='lower')
        ax[index].set_xlabel('Longitude')
        ax[index].set_ylabel('Latitude')
        ax[index].set_title(r'$lambda_0$={}'.format(lamda0))
        ax[index].set_yticks([0,5,10,15,20])
        ax[index].set_yticklabels(['25','30','35','40','45','50'])
        ax[index].set_xticks([2,12,22,32,42,52])
        ax[index].set_xticklabels(['-125','-115','-105','-95','-85','-75'])

    # plot misfitreduce as a function of lamda0
    ax[-1].scatter(lamda0_list,misfitreduce,s=3,color='k')
    ax[-1].set_xlabel(r'$lambda_0$')
    ax[-1].set_ylabel('Reduced misfit (%)')
    fig.tight_layout()
    plt.savefig('{}/solve_matrix_onedep_lamda.pdf'.format(savepath))
    plt.colorbar(cm,ax=ax)
    plt.close()

def main():
    filename = 'rp_SAPA_amp_400.BHT.vel.snrS.xh'
    dirname = '/home/meichen/work1/SH_TOMO/Combined_xh'
    curpath = '/home/meichen/Research/SH_TOMO/USArray'

    readbinary(filename=filename,dirname=dirname,savepath=curpath)

main()
