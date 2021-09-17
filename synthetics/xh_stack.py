#!/home/meichen/anaconda3/bin/python
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
        dt = kwargs.get('dt')
        dist_min = kwargs.get('dist_min')
        dist_max = kwargs.get('dist_max')
        #rtfilename = kwargs.get('rtfilename')

    stack_seism = np.zeros(int((cut_b+cut_a)/dt))
    stack_num = 0
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

                if plot == 'True' and flt[6]>dist_min and flt[6]<dist_max:

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
                    N2 = N1 + int((cut_a+cut_b)/delta)
                    t = np.linspace(-1*cut_b,cut_a,num=len(amp[N1:N2]))
                    if normalize == 'True':
                        norm = np.array(amp[N1:N2])/amp[N0]*0.1
                        ax.plot(t,norm+yy_value,c='k',lw=0.5,alpha=0.5)
                        norm = norm / np.sqrt(np.abs(norm)) # root stack
                        stack_seism = stack_seism + norm
                        stack_num = stack_num + 1
            else:
                break
    fp.close()

    stack_seism= np.array(stack_seism)
    stack_seism = stack_seism / np.float(stack_num)
    stack_seism = root_square(stack_seism)

    # plot stack
    ax.plot(t,stack_seism+dist_min-0.1,c='r',lw=1,alpha=0.8)
    ax.tick_params('x',bottom=True,top=True)
    ax.set_xticks([-50,0,50,100,150,200,250,300])
    ax.set_xticklabels(['-50','0','50','100','150','200','250','300'])
    #ax.set_ylim([ylim_max,ylim_min])
#    ax.set_xlim([-50,300])
    ax.set_xlabel('Time (s) aligned to {}'.format(phase))
    ax.set_ylabel('{}'.format(ylabel))
    ax.set_title('distance deg {}-{}'.format(dist_min,dist_max))

def main():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    path = '/home/meichen/work1/SH_TOMO/synthetics/XH_files/PREM_NOCRUST_10s'
    current_path = '/home/meichen/Research/SH_TOMO/synthetics'
    
    dirname = '{}'.format(path)
    fig,ax = plt.subplots(1,1,figsize=[7,4])
    filename = "prem_nocrust_10s.sflip.BXT.dis.xh"
    readxh(ax=ax,filename=filename,dirname=dirname,odirname=current_path,plot='True',align_index=3,phase='S',yy_index=6,ylabel='Distance (deg)',cut_b=50,cut_a=300,ofile='{}.pdf'.format(filename),normalize='True',dt=0.063,dist_min=80.0,dist_max=80.1)
    plt.savefig('{}/{}.diststack.pdf'.format(current_path,filename))
    plt.close()

main()
