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

                ax1.fill_between(t,np.array(seism)*0.2+np.float(gcarc),np.float(gcarc),where=(np.array(seism) > 0), facecolor='b',alpha=0.5)
                ax1.fill_between(t,np.array(seism)*0.2+np.float(gcarc),np.float(gcarc),where=(np.array(seism) < 0), facecolor='r',alpha=0.5)
#                ax1.plot(t,np.array(seism)+np.float(gcarc),lw=0.5,ls='-',alpha=0.8,color='k')

#                ax1.plot(0,seism[3499]+gcarc,'ro',markersize=2)
#                ax1.plot(t[hintg1],seism[hintg1]+rp_dist,'bo',markersize=2)
#                ax1.plot(t[hintg2],seism[hintg2]+rp_dist,'mo',markersize=2)
#                ax1.plot(t[hintg4],seism[hintg4]+rp_dist,'go',markersize=2)
#                ax1.plot(t[hintg5],seism[hintg5]+rp_dist,'yo',markersize=2)
                seism_num = seism_num + 1
                seism_root = root_stack(seism,2)
                print(t_to_S,np.max(seism),np.max(seism_root))
                stack = stack + np.array(seism_root)
            else:
                break

    ax1.set_xlabel('Time (s)',size=26)
    ax1.set_ylabel('Distance (deg)',size=26)
    ax1.tick_params('both',labelsize=26)
    ax1.set_title('{} num={}'.format(filename,seism_num),size=30)
    ax1.set_xlim([-250,250])
    stack = stack / seism_num
    env = obspy.signal.filter.envelope(stack)
#    stack = stack / np.amax(env)
#    env = env/np.amax(env)
#    ax2.fill_between(np.arange(7000)*0.1-350,stack,0,where=(stack>0),facecolor='b')
#    ax2.fill_between(np.arange(7000)*0.1-350,stack,0,where=(stack<0),facecolor='r')
    ax2.plot(np.arange(7000)*0.1-350,stack,'k-',lw=0.5)
#    ax2.plot(np.arange(7000)*0.1-350,env,color='k',linestyle='-',lw=1)
#    ax2.plot(np.arange(100)*0.1-5,env[3450:3550],'b-',lw=0.5)
    ax2.tick_params('both',labelsize=26)
    ax2.set_xlabel('Time (s)',size=26)
    ax2.set_ylabel('Amplitude',size=26)
    ax2.set_title('Stacked',size=30)
    ax2.set_xlim([-250,250])
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
