#!/home/meichen/anaconda3/bin/python

def root_stack(seism):                                                      
    import numpy as np
    seism_root = np.zeros(len(seism))
    for i in np.arange(len(seism)):
        sign = seism[i]/np.abs(seism[i])
        seism_root[i] = sign*np.sqrt(np.abs(seism[i]))
    return seism_root

def root_square(seism):
    import numpy as np
    seism_square = np.zeros(len(seism))
    for i in np.arange(len(seism)):
        sign = seism[i]/np.abs(seism[i])
        seism_square[i] = sign*seism[i]*seism[i]
    return seism_square

def readxh(**kwargs):

##---------------------------------------##

# This function read the xh file header.

# Created by Meichen Liu on July 6th, 2019.

##---------------------------------------##
##parameters
# filename	the name of the xh file
# dirname	the directory of filename
# plot		True or False
##if plot is "True", following parameters are needed
# align_index	tpck[0-19]
# yy_index	[0-19]
# cut_b		the window to plot start before tpck (sec)
# cut_a		the window to plot start after tpck (sec)
# odirname	the directory to save output file
# ofile		the name of output file
# normalize	True or False
# phase		the name of aligned phase
# ylabel	the name of ylabel
# rd		reflection depth (km) of the reverberation phase


    import struct
    from obspy.core import UTCDateTime
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    import matplotlib.cm as cms
    import matplotlib.colors as colors

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    plot = kwargs.get('plot')
    if plot == 'True':
        ax = kwargs.get('ax')
        align_index = kwargs.get('align_index')
        yy_index = kwargs.get('yy_index')
        cut_b = kwargs.get('cut_b')
        cut_a = kwargs.get('cut_a')
        ofile = kwargs.get('ofile')
        odirname = kwargs.get('odirname')
        normalize = kwargs.get('normalize')
        phase = kwargs.get('phase')
        ylabel = kwargs.get('ylabel')
        color_seism = kwargs.get('color_seism')
        rd = kwargs.get('rd')
        rtfilename = kwargs.get('rtfilename')
    
    os.chdir('{}'.format(dirname))
    sw = np.zeros(int((cut_b+cut_a)/0.1)-1)
    se = np.zeros(int((cut_b+cut_a)/0.1)-1)
    nw = 0.0
    ne = 0.0

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
                evtcd = head[217]
                netw = head[218]
                stnm = head[219]
                chan = head[220]
                rcomment = head[221]
                wavf = head[222]
                padding = head[223]
        
                xhdata = fp.read(ndata*4)
                xhdata_format = "f" * ndata
                amp = struct.unpack(xhdata_format,xhdata)

                if plot == 'True' and slat>36 and slat<41 and slon>-120 and slon<-105:

                    print(evtcd)
                    ot_utc = UTCDateTime("{}-{}-{}T{}:{}:{}".format(ot[0],ot[1],ot[2],ot[3],ot[4],ot[5]))
                    tstart_utc = UTCDateTime("{}-{}-{}T{}:{}:{:.2f}".format(tstart[0],tstart[1],tstart[2],tstart[3],tstart[4],tstart[5]))
                    b = tstart_utc - ot_utc
                    e = b + delta*ndata
                    align_time = tpck[int(align_index)]

                    yy_value = flt[yy_index]

                    if b > align_time - cut_b or e < align_time + cut_a:
                        print("Time window exceeds the range of trace {}.{} ...".format(netw,stnm))
                        continue
                    N0 = int((-1*b + align_time)/delta)
                    N1 = int((-1*b + align_time - cut_b)/delta)
                    N2 = int((-1*b + align_time + cut_a)/delta)
                    t = np.linspace(-1*cut_b,cut_a,num=len(amp[N1:N2]))
                    if normalize == 'True' and slon<-100:
                        norm = np.array(amp[N1:N2])/amp[int((-b+tpck[3])/delta)]
#                        norm = norm/np.sqrt(np.abs(norm))  # root stack
                        sw = sw + norm[0:len(sw)]
                        nw = nw + 1
                        ax.plot(t,norm+yy_value,lw=0.1,alpha=0.1,c='k')
                    elif normalize == 'True' and slon>=-100:
                        norm = np.array(amp[N1:N2])/amp[int((-b+tpck[3])/delta)]
#                        norm = norm/np.sqrt(np.abs(norm))  # root stack
                        se = se + norm[0:len(se)]
                        ne = ne + 1
                        ax.plot(t,norm+yy_value,lw=0.1,alpha=0.1,c='k')
                    elif normalize == 'False':
                        ax.plot(t,np.array(amp[N1:N2])+yy_value,lw=0.1,alpha=0.4)
            else:
                break
    fp.close()

    print(nw)
    print(ne)
    sw = sw / nw
    se = se / ne
    fig1,ax1 = plt.subplots(1,1,figsize=[7,3])
    ax1.plot(t[0:len(sw)],sw,lw=1,c='k',label='west')
    ax1.plot(t[0:len(se)],se,lw=1,c='r',label='east')
    ax1.set_ylim([-0.2,0.2])
    ax1.set_xlim([-50,300])
    ax1.set_xlabel('Time (s)',size=12)
    ax1.set_ylabel('Amp',size=12)
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig('eastwest.pdf')

    ax.set_ylim([110,60])
#    ax.set_xlim([-50,400])
    ax.set_xlabel('Time (s) aligned to {}'.format(phase))
    ax.set_ylabel('{}'.format(ylabel))
#    # plot 410 and 660 reflection phase curves
#    depth1 = 400
#    depth2 = 670
#    depth3 = 250
#    rt = np.genfromtxt('{}'.format(rtfilename))
#    reflection_depth = rt[:,0]
#    event_depth = rt[:,1]
#    distance = rt[:,2]
#    time_resp_S = rt[:,4]
#
#    event_depth_410 = event_depth[np.abs(reflection_depth-depth1)<0.1]
#    time_resp_S_410 = time_resp_S[np.abs(reflection_depth-depth1)<0.1]
#    distance_410 = distance[np.abs(reflection_depth-depth1)<0.1]
#    time_resp_S_410 = time_resp_S_410[np.abs(event_depth_410-20)<0.1]
#    distance_410 = distance_410[np.abs(event_depth_410-20)<0.1]
#    distance_410 = distance_410[time_resp_S_410>0.0]
#    time_resp_S_410 = time_resp_S_410[time_resp_S_410>0.0]
#    ax.plot(time_resp_S_410*10+500,distance_410*10,color='k',lw=0.1,alpha=0.5)
#    event_depth_660 = event_depth[np.abs(reflection_depth-depth2)<0.1]
#    time_resp_S_660 = time_resp_S[np.abs(reflection_depth-depth2)<0.1]
#    distance_660 = distance[np.abs(reflection_depth-depth2)<0.1]
#    time_resp_S_660 = time_resp_S_660[np.abs(event_depth_660-20)<0.1]
#    distance_660 = distance_660[np.abs(event_depth_660-20)<0.1]
#    distance_660 = distance_660[time_resp_S_660>0.0]
#    time_resp_S_660 = time_resp_S_660[time_resp_S_660>0.0]
#    ax.plot(time_resp_S_660*10+500,distance_660*10,color='k',lw=0.1,alpha=0.5)
#    event_depth_250 = event_depth[np.abs(reflection_depth-depth3)<0.1]
#    time_resp_S_250 = time_resp_S[np.abs(reflection_depth-depth3)<0.1]
#    distance_250 = distance[np.abs(reflection_depth-depth3)<0.1]
#    time_resp_S_250 = time_resp_S_250[np.abs(event_depth_250-20)<0.1]
#    distance_250 = distance_250[np.abs(event_depth_250-20)<0.1]
#    distance_250 = distance_250[time_resp_S_250>0.0]
#    time_resp_S_250 = time_resp_S_250[time_resp_S_250>0.0]
#    ax.plot(time_resp_S_250*10+500,distance_250*10,color='k',lw=0.1,alpha=0.5)
    
#    plt.colorbar(cm,ax=ax)

#    deg = table[:,0]
#    tS = table[:,1]
#    tScS = table[:,2]
#    tSSv410s = table[:,3]
#    tsS = table[:,4]
#    tSSv660s = table[:,5]
#    tsScS = table[:,6]
#    tSS = table[:,7]
#    tS410S = table[:,8]
#    tS660S = table[:,9]
#    ax.plot(tScS-tS,deg,color='k',lw=0.2)
#    ax.plot(tSSv410s-tS,deg,color='k',lw=0.2)
#    ax.plot(tsS-tS,deg,color='k',lw=0.2)
#    ax.plot(tSSv660s-tS,deg,color='k',lw=0.2)
#    ax.plot(tsScS-tS,deg,color='k',lw=0.2)
#    ax.plot(tSS-tS,deg,color='k',lw=0.2)
#    ax.plot(tS410S-tS,deg,color='k',lw=0.2)
#    ax.plot(tS660S-tS,deg,color='k',lw=0.2)

def main():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    path = '/home/meichen/work1/SH_TOMO/Combined_xh/finalxhfiles'
    current_path = '/home/meichen/Research/SH_TOMO/USArray'
    
    dirname = '{}'.format(path)
    fig,ax = plt.subplots(1,1,figsize=[7,11])
    filename = "all.BHT.dis.rmvtrace.xh.gcarc60110"
    readxh(ax=ax,filename=filename,dirname=dirname,odirname=current_path,plot='True',align_index=3,phase='S',yy_index=6,ylabel='Distance (deg)',cut_b=50,cut_a=300,ofile='{}.pdf'.format(filename),normalize='True',color_seism='black')
    fig.savefig('{}/{}.wf.pdf'.format(current_path,filename))
    plt.close()

main()
