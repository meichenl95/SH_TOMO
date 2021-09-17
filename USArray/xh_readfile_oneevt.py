#!/home/meichen/anaconda3/bin/python

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
# yy		flt or intg
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

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    plot = kwargs.get('plot')
    if plot == 'True':
        ax = kwargs.get('ax')
        align_index = kwargs.get('align_index')
        yy = kwargs.get('yy')
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
        evtcode = kwargs.get('evtcd')

    os.chdir('{}'.format(dirname))

    xhhead_format = "fii" + "f"*11 + "iiiiif"*2 + "ifffiii" + "ff"*60 + "fff"+"f"*20*2+"i"*20 + "14s8s8s8s8s72s8s34s"
    root_stack = np.zeros((int) ((cut_a+cut_b)/0.1))
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
                evtcd = ''.join(list(head[217].decode())[0:7])
                netw = head[218].decode()
                stnm = head[219].decode()
                chan = head[220]
                rcomment = head[221]
                wavf = head[222]
                padding = head[223]
        
                xhdata = fp.read(ndata*4)
                xhdata_format = "f" * ndata
                amp = struct.unpack(xhdata_format,xhdata)

                if plot == 'True' and evtcd == evtcode:

                    print('yes')
                    ot_utc = UTCDateTime("{}-{}-{}T{}:{}:{}".format(ot[0],ot[1],ot[2],ot[3],ot[4],ot[5]))
                    tstart_utc = UTCDateTime("{}-{}-{}T{}:{}:{}".format(tstart[0],tstart[1],tstart[2],tstart[3],tstart[4],tstart[5]))
                    b = tstart_utc - ot_utc
                    e = b + delta*ndata
                    align_time = tpck[int(align_index)]
                    yy_value = flt[yy_index]
                    

                    if b > align_time - cut_b or e < align_time + cut_a:
                        print("Time window exceeds the range of trace {}.{} in event {} ...".format(netw,stnm,evtcd))
                        continue
                    N0 = int((-1*b + align_time)/delta)
                    N1 = int((-1*b + align_time - cut_b)/delta)
                    N2 = int((cut_b + cut_a)/0.1) + N1
                    t = np.linspace(-1*cut_b,cut_a,num=len(amp[N1:N2]))
                    if normalize == 'True':
                        if amp[N0] - min(amp[N1:N2]) != 0:
                            norm = np.array(amp[N1:N2])/amp[N0]
#                            norm = norm/np.sqrt(np.abs(norm))
                            ax.plot(t,norm + yy_value,lw=0.4,c='{}'.format(color_seism),alpha=0.5)
#                            ax.scatter(t[intg[1]-N1],amp[intg[1]]/(max(amp[N1:N2])-min(amp[N1:N2]))*2 + yy_value,marker='o',color='blue',s=1,alpha=0.6)
#                            ax.scatter(t[intg[2]-N1],amp[intg[2]]/(max(amp[N1:N2])-min(amp[N1:N2]))*2 + yy_value,marker='o',color='green',s=1,alpha=0.6)
#                            ax.scatter(t[intg[4]-N1],amp[intg[4]]/(max(amp[N1:N2])-min(amp[N1:N2]))*2 + yy_value,marker='o',color='yellow',s=1,alpha=0.6)
#                            ax.scatter(t[intg[5]-N1],amp[intg[5]]/(max(amp[N1:N2])-min(amp[N1:N2]))*2 + yy_value,marker='o',color='orange',s=1,alpha=0.6)
                    elif normalize == 'False':
                        ax.plot(t,np.array(amp[N1:N2])/amp[(int) ((tpck[3]-b)/delta)] + yy_value,lw=0.4,c='{}'.format(color_seism),alpha=0.3)
                        #ax.plot(t,np.array(amp[N1:N2])/np.sqrt(np.abs(amp[N1:N2])) + yy_value,lw=0.4,c='{}'.format(color_seism),alpha=0.3)
                        root_stack = root_stack + np.array(amp[N1:N2])/np.sqrt(np.abs(amp[N1:N2]))
#                        ax.scatter(0,amp[int((-1*b+align_time)/delta)] + yy_value,marker='o',color='red',s=3)
                        ax.scatter(tpck[3]-align_time,amp[int((-1*b+tpck[3])/delta)] + yy_value,marker='o',color='green',s=3)
                    else:
                        print("Please specify normalize option: True or False?")
            else:
                break
    ax.tick_params('x',bottom=True,top=True)
    ax.set_xlim([-1*cut_b,cut_a])
    ax.set_xlabel('Time (s) aligned to {}'.format(phase))
    ax.set_ylabel('{}'.format(ylabel))
    fp.close()

def single_event():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    path = '/home/meichen/work1/SH_TOMO/Combined'
    current_path = '/home/meichen/Research/SH_TOMO/USArray'
    
    events_selected = pd.read_csv('{}/events_selected.txt'.format(current_path),skipinitialspace=True,header=None,sep=' ')
    events_selected = np.array(events_selected)
    dirname = '/home/meichen/work1/SH_TOMO/Combined_xh'
    fig,ax = plt.subplots(1,1,figsize=[9,4])
    filename = "all_MidAtl.BHT.vel.lowsnrS.xh"
    evtcode = '120909D'
    readxh(ax=ax,filename=filename,dirname=dirname,odirname=current_path,plot='True',align_index=3,phase='S',yy='intg',yy_index=6,ylabel='Distance (deg)',cut_b=300,cut_a=100,ofile='{}.pdf'.format(filename),normalize='True',color_seism='black',evtcd=evtcode)
    fig.savefig('{}/{}.pdf'.format(current_path,evtcode))
    plt.close()

single_event()
