#!/home/meichen/anaconda3/bin/python

def root_stack(seism,n):

    import numpy as np

    seism_root = np.zeros(len(seism))
    for i in np.arange(len(seism)):
        sign = seism[i]/np.abs(seism[i])
        seism_root[i] = sign*np.sqrt(np.abs(seism[i]))
    return seism_root

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
    import obspy.signal.filter

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    savepath = kwargs.get('savepath')
    ndata = kwargs.get('ndata')
    
    seism_num = 0
    stack_seism = np.zeros((1800,7000))
    stack_num = np.zeros(1800)
    stack = np.zeros(7000)

    os.chdir('{}'.format(dirname))
    
    fig = plt.figure(constrained_layout=True,figsize=[14,24])
    gs = fig.add_gridspec(6,1)
    ax1 = fig.add_subplot(gs[0:5,0])
    ax2 = fig.add_subplot(gs[5,0])
    bihead_format = "f"*7003 + "i"*4
    with open("{}".format(filename),"rb") as fp:
        while True:
            bihead = fp.read(28028)
            if bihead:
                head = struct.unpack(bihead_format,bihead)

                seism = head[0:7000]
                t_to_S = head[7000]
                t = np.arange(7000)*0.1 -350
                rp_dist = head[7001]
                gcarc = head[7002]
                hintg1 = head[7003] + 3499
                hintg2 = head[7004] + 3499
                hintg4 = head[7005] + 3499
                hintg5 = head[7006] + 3499
                seism_root = root_stack(seism,2)

                stack_seism[np.int(gcarc*10),:] = stack_seism[np.int(gcarc*10),:] + seism_root
                stack_num[np.int(gcarc*10)] = stack_num[np.int(gcarc*10)] + 1
                seism_num = seism_num + 1
                print(t_to_S,np.max(seism),np.max(seism_root),np.int(gcarc*10))
                stack = stack + np.array(seism_root)
            else:
                break

    ylim_min = np.where(stack_num>0)[0][0]
    ylim_max = np.where(stack_num>0)[0][-1]
    for i in np.arange(1800):
        if stack_num[i]>0:
            stack_seism[i,:] = stack_seism[i,:]/np.float(stack_num[i])

    cm = ax1.imshow(stack_seism,cmap='bwr',vmin=-1.0,vmax=1.0,aspect='auto')
    ax1.tick_params('x',bottom=True,top=True)
    ax1.set_xlabel('Time (s)',size=26)
    ax1.set_ylabel('Distance (deg)',size=26)
    ax1.tick_params('both',labelsize=26)
    ax1.set_title('{} num={}'.format(filename,seism_num),size=30)
    ax1.set_xticks([1500,2500,3500,4500,5500,6500])
    ax1.set_xticklabels(['-200','-100','0','100','200','300'])
    ax1.set_ylim([ylim_min,ylim_max])
    labels = ax1.get_yticks().tolist()
    labels = [np.float(labels[i])/10.0 for i in np.arange(len(labels))]
    ax1.set_yticklabels(list(labels))
    stack = stack / seism_num
    ax2.plot(np.arange(7000)*0.1-350,stack,'k-',lw=0.5)
    ax2.tick_params('both',labelsize=26)
    ax2.set_xlabel('Time (s)',size=26)
    ax2.set_ylabel('Amplitude',size=26)
    ax2.set_title('Stacked',size=30)
    ax2.set_xlim([-250,250])
    ax1.set_xlim([1000,6000])
    ax2.hlines(0,xmin=-250,xmax=250,color='r',ls='--')
    plt.savefig('{}/{}.png'.format(savepath,filename))
    print(np.max(stack))
    max_seism = 0.0
    for i in np.arange(3475,3525):
        if stack[i]**2>max_seism**2:
            max_i = i
    print(stack[i])


def main():
    filename = 'seism_400_800.BHT.vel.peakS.xh.5.780.32.-120.5'
    dirname = '/home/meichen/work1/SH_TOMO/Combined_xh'
    curpath = '/home/meichen/Research/SH_TOMO/USArray'

    readbinary(filename=filename,dirname=dirname,savepath=curpath)

main()
