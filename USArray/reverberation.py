#!/home/meichen/anaconda3/bin/python


def readbinary(**kwargs):

# parameters
# filename		the name of the binary file to read in
# dirname		the directory of the binary file
# savepath		the directory to save output
# ndata			the number of data in each seismogra
# stack_num		the minimum number of stacked seismograms of each grid point

    import struct
    import numpy as np
    import os
    import numpy as np
    import matplotlib.cm as cmx
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from mpl_toolkits.mplot3d import Axes3D

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    savepath = kwargs.get('savepath')
    ndata = kwargs.get('ndata')
    stack_num = kwargs.get('stack_num')
    

    os.chdir('{}'.format(dirname))
    depth = []
    lon = []
    lat = []
    amp = []
    seism_num = []
    
    bihead_format = "f"*5
    with open("{}".format(filename),"rb") as fp:
        while True:
            bihead = fp.read(20)
            if bihead:
                head = struct.unpack(bihead_format,bihead)

                depth.append(head[0])
                lon.append(head[1])
                lat.append(head[2])
                seism_num.append(head[4])
                if(head[4]>stack_num):
                    amp.append(head[3])
                    print(head[0],head[1],head[2],head[3])
                else:
                    amp.append(0)
            else:
                break

    depth = np.array(depth)*(-1.0)
    lon = np.array(lon)
    lat = np.array(lat)
    amp = np.array(amp)
    seism_num = np.array(seism_num)

    # 3d figure
    fig = plt.figure(figsize=[8,6])
    ax = fig.gca(projection='3d')
    
    norm = colors.DivergingNorm(vmin=-0.5,vmax=0.5,vcenter=0)
    cmap = cmx.ScalarMappable(norm=norm,cmap=plt.get_cmap('seismic'))
    cmap.set_array([])
    ax.scatter(lon,lat,depth,c=cmap.to_rgba(amp),marker='o',s=1,alpha=0.5)
    plt.colorbar(cmap,ax=ax)
    ax.set_xlabel('longitude',size=14)
    ax.set_ylabel('latitude',size=14)
    ax.set_zlabel('depth',size=14)
    
    fig.savefig('{}/reverberation.{}.pdf'.format(savepath,filename))

    # amp figure
    fig,ax = plt.subplots(2,10,figsize=[70,6])
    mrshp = amp.reshape(80,61,26)
    for i in np.arange(10):
        cm = ax[0][i].imshow(mrshp[:,i*5+7,:],cmap='seismic_r',vmin=-0.5,vmax=0.5,aspect='auto')
        ax[0][i].set_title('lon={}'.format(i*5+7-127))
        ax[0][i].set_xlabel('latitude',size=10)
        ax[0][i].set_xticks([3,5,7,9,11,13,15,17,19,21,23,25,27])
        ax[0][i].set_yticks([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80])
        ax[0][i].set_xticklabels(['','30','','','','','40','','','','','50',''])
        ax[0][i].set_yticklabels(['','','','200','','','','400','','','','600','','','','800'])
        ax[0][i].tick_params(axis='both',labelsize=10)
        ax[0][i].set_xlim([2,27])
        plt.colorbar(cm,ax=ax[0][i])
        cm = ax[1][i].imshow(mrshp[:,:,i*2+1],cmap='seismic_r',vmin=-0.5,vmax=0.5,aspect='auto')
        ax[1][i].set_title('lat={}'.format(i*2+1+25))
        ax[1][i].set_xlabel('logitude',size=10)
        ax[1][i].set_xticks([-3,2,7,12,17,22,27,32,37,42,47])
        ax[1][i].set_yticks([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80])
        ax[1][i].set_xticklabels(['','','-120','','','','-100','','','','-80'])
        ax[1][i].set_yticklabels(['','','','200','','','','400','','','','600','','','','800'])
        ax[1][i].tick_params(axis='both',labelsize=10)
        ax[1][i].set_xlim([-3,47])
        plt.colorbar(cm,ax=ax[1][i])
    fig.savefig('{}/profiles.{}.pdf'.format(savepath,filename))

    # count figure
    fig,ax = plt.subplots(2,10,figsize=[70,6])
    mrshp = seism_num.reshape(80,61,26)
    for i in np.arange(10):
        cm = ax[0][i].imshow(mrshp[:,i*5+7,:],cmap='rainbow',vmin=0,vmax=np.max(seism_num),aspect='auto')
        ax[0][i].set_title('lon={}'.format(i*5+7-127))
        ax[0][i].set_xlabel('latitude',size=10)
        ax[0][i].set_xticks([3,5,7,9,11,13,15,17,19,21,23,25,27])
        ax[0][i].set_yticks([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80])
        ax[0][i].set_xticklabels(['','30','','','','','40','','','','','50',''])
        ax[0][i].set_yticklabels(['','','','200','','','','400','','','','600','','','','800'])
        ax[0][i].tick_params(axis='both',labelsize=10)
        ax[0][i].set_xlim([2,27])
        plt.colorbar(cm,ax=ax[0][i])
        ax[1][i].imshow(mrshp[:,:,i*2+1],cmap='rainbow',vmin=0,vmax=np.max(seism_num),aspect='auto')
        ax[1][i].set_title('lat={}'.format(i*2+1+25))
        ax[1][i].set_xlabel('logitude',size=10)
        ax[1][i].set_xticks([-3,2,7,12,17,22,27,32,37,42,47])
        ax[1][i].set_yticks([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80])
        ax[1][i].set_xticklabels(['','','-120','','','','-100','','','','-80'])
        ax[1][i].set_yticklabels(['','','','200','','','','400','','','','600','','','','800'])
        ax[1][i].tick_params(axis='both',labelsize=10)
        ax[1][i].set_xlim([-3,47])
    fig.savefig('{}/counts.{}.pdf'.format(savepath,filename))
    plt.close()

def main():
    filename = 'rp_400_800.BHT.vel.peakS.xh.5d.5s'
    dirname = '/home/meichen/work1/SH_TOMO/Combined_xh'
    curpath = '/home/meichen/Research/SH_TOMO/USArray'

    readbinary(filename=filename,dirname=dirname,savepath=curpath,stack_num=100)

main()
