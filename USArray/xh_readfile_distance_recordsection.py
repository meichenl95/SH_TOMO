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
        seism_square[i] = sign*np.sqrt(abs(seism[i]))
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
        dt = kwargs.get('dt')

    stack_seism = np.zeros((1800,int((cut_a+cut_b)/dt)))
    stack_num = np.zeros(1800)
    # mark reverberation phases
    
    os.chdir('{}'.format(dirname))

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

                if plot == 'True' :

                    ot_utc = UTCDateTime("{}-{}-{}T{}:{}:{}".format(ot[0],ot[1],ot[2],ot[3],ot[4],ot[5]))
                    tstart_utc = UTCDateTime("{}-{}-{}T{}:{}:{:.2f}".format(tstart[0],tstart[1],tstart[2],tstart[3],tstart[4],tstart[5]))
                    b = tstart_utc - ot_utc
                    e = b + delta*ndata
                    align_time = tpck[int(align_index)]

                    yy_value = flt[yy_index]*10

                    if b > align_time - cut_b or e < align_time + cut_a:
                        print("Time window exceeds the range of trace {}.{} of event {} ...".format(netw,stnm,evtcd))

                        continue
                    N0 = int((-1*b + align_time)/delta)
                    N1 = int((-1*b + align_time - cut_b)/delta)
                    N2 = N1 + int((cut_b+cut_a)/delta)+1
                    t = np.linspace(-1*cut_b,cut_a,num=len(amp[N1:N2]))
                    if normalize == 'True':
                        norm = np.array(amp[N1:N2])/amp[N0]
#                        norm = root_square(norm)  # root stack
                        stack_seism[np.int(yy_value),:] = stack_seism[np.int(yy_value),:] + norm
                        stack_num[np.int(yy_value)] = stack_num[np.int(yy_value)] + 1
            else:
                break
    fp.close()

    print(np.sum(stack_num))
    ylim_min = np.where(stack_num>0)[0][0]
    ylim_max = np.where(stack_num>0)[0][-1]
    for i in np.arange(1800):
        if stack_num[i] > 0:
            stack_seism[i,:] = stack_seism[i,:] / np.float(stack_num[i])
#            stack_seism[i,:] = stack_seism[i,:] * stack_seism[i,:]
#            stack_seism[i,:] = root_square(stack_seism[i,:])
#            for j in np.arange(3500):
#            ax.fill_between(np.arange(3500)*0.1-50,stack_seism[i,:]+i*0.1+70, i*0.1+70,where=(stack_seism[i,:] > 0), facecolor='r')
#            ax.fill_between(np.arange(3500)*0.1-50, stack_seism[i,:]+i*0.1+70, i*0.1+70,where=(stack_seism[i,:] < 0), facecolor='b')

    # plot record section
    cm = ax.imshow(stack_seism,cmap='bwr',vmin=-0.1,vmax=0.1,aspect='auto')
    ax.tick_params('x',bottom=True,top=True)
    ax.set_xticks([0,500,1000,1500,2000,2500,3000,3500,4000])
    ax.set_xticklabels(['-50','0','50','100','150','200','250','300','350'])
    ax.set_yticks([0,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700])
    ax.set_yticklabels(['0','30','40','50','60','70','80','90','100','110','120','130','140','150','160','170'])
    #ax.set_ylim([ylim_max,ylim_min])
    ax.set_ylim([1100,600])
#    ax.set_xlim([-50,300])
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
    
    plt.colorbar(cm,ax=ax)

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
    filename = "all.BHT.dis.snrS2.xh.gcarc60110.20evt"
    rtfilename = "Utils/TauP_calculations/TravelTimeTables/rt_dep5.prem_ml..s"
    readxh(ax=ax,filename=filename,dirname=dirname,odirname=current_path,plot='True',rtfilename=rtfilename,align_index=3,phase='S',yy_index=6,ylabel='Distance (deg)',cut_b=50,cut_a=400,ofile='{}.pdf'.format(filename),normalize='True',color_seism='black',dt=0.1)
    plt.savefig('{}/{}.pdf'.format(current_path,filename))
    plt.close()

main()
