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

    # mark reverberation phases
    os.chdir('/home/meichen/Utils/TauP_calculations/TravelTimeTables/')
    rt = np.genfromtxt('rt.prem_ml..s')
    
    os.chdir('{}'.format(dirname))

    xhhead_format = "fii" + "f"*11 + "iiiiif"*2 + "ifffiii" + "ff"*60 + "fff"+"f"*20*2+"i"*20 + "14s8s8s8s8s72s8s34s"
    num = 0
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

                if plot == 'True':

                    ot_utc = UTCDateTime("{}-{}-{}T{}:{}:{}".format(ot[0],ot[1],ot[2],ot[3],ot[4],ot[5]))
                    tstart_utc = UTCDateTime("{}-{}-{}T{}:{}:{}".format(tstart[0],tstart[1],tstart[2],tstart[3],tstart[4],tstart[5]))
                    b = tstart_utc - ot_utc
                    e = b + delta*ndata
                    align_time = tpck[int(align_index)]
                    yy_value = flt[yy_index]
                    

                    if b > align_time - cut_b or e < align_time + cut_a:
                        print("Time window exceeds the range of trace {}.{} ...".format(netw,stnm))
                        continue
                    N0 = int((-1*b + align_time)/delta)
                    N1 = int((-1*b + align_time - cut_b)/delta)
                    N2 = int((cut_b + cut_a)/delta) + N1
                    t = np.arange(N1-N0,N2-N0)*delta
                    if normalize == 'True':
                        norm = np.array(amp[N1:N2])/amp[N0]
#                        norm = norm/np.sqrt(np.abs(norm))
                        ax.plot(t,norm + yy_value,lw=0.4,c='{}'.format(color_seism),alpha=0.5)
                    elif normalize == 'False':
                        ax.plot(t,np.array(amp[N1:N2]) ,lw=0.4,c='{}'.format(color_seism),alpha=0.3)
                    else:
                        print("Please specify normalize option: True or False?")
                    if (np.abs(num)<0.001):
                        stack = np.array(amp[N1:N2])
                    else:
                        stack = stack + np.array(amp[N1:N2])
                    num = num + 1
            else:
                break
    ax.plot(t,stack/num,lw=0.5,c='r')
    ax.tick_params('x',bottom=True,top=True)
    ax.set_xlim([-1*cut_b,cut_a])
    #ax.set_ylim([50,110])
    ax.set_xlabel('Time (s) aligned to {}'.format(phase))
    ax.set_ylabel('{}'.format(ylabel))
    fp.close()

def main():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    path = '/home/meichen/work1/SH_TOMO/synthetics/XH_files/MTZ_2X2_NOCRUST_10s'
    current_path = '/home/meichen/Research/SH_TOMO/synthetics'
    
    events_selected = pd.read_csv('{}/events_cat.txt'.format(current_path),skipinitialspace=True,header=None,sep=' ')
    events_selected = np.array(events_selected)
#    for eventid in events_selected[0:1,0]:
    for eventid in ['010100X']:
        dirname = '{}'.format(path)
        fig,ax = plt.subplots(1,1,figsize=[5,12])
        filename = "{}.dis.xh".format(eventid)
        readxh(ax=ax,filename=filename,dirname=dirname,odirname=current_path,plot='True',align_index=1,phase='S',yy='intg',yy_index=6,ylabel='Distance (deg)',cut_b=50,cut_a=50,ofile='{}.pdf'.format(filename),normalize='False',color_seism='black')
        plt.savefig('{}/{}.pdf'.format(current_path,filename))
        plt.close()

def single_gridpoint():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    current_path = '/home/meichen/Research/SH_TOMO/synthetics'
    
    dirname = '/home/meichen/work1/SH_TOMO/synthetics/XH_files/pointcheck'
    fig,ax = plt.subplots(1,1,figsize=[9,4])
    filename = "prem_nocrust_10s_rsmp.bp_102070100.snrS.BXT.vel.xh.dep_665.deg_1.lat35.lon-115.win5"
    readxh(ax=ax,filename=filename,dirname=dirname,odirname=current_path,plot='True',align_index=9,phase='SSvDEPs',yy='flt',yy_index=9,ylabel='Distance (deg)',cut_b=2.5,cut_a=2.5,ofile='{}.pdf'.format(filename),normalize='False',color_seism='black')
    fig.savefig('{}/{}.pdf'.format(current_path,filename))
    plt.close()

def plot_overlap():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    path = '/home/meichen/work1/SH_TOMO/events'
    current_path = '/home/meichen/Research/SH_TOMO/USArray'
    
    events_selected = pd.read_csv('{}/events_selected.txt'.format(current_path),skipinitialspace=True,header=None,sep=' ')
    events_selected = np.array(events_selected)
    for eventid in events_selected[:,0]:
        dirname = '{}/event_{}/waveforms/XH_files'.format(path,eventid)
        fig,ax = plt.subplots(1,1,figsize=[5,9])
        filename1 = "{}.BHT.dis.snrS.xh".format(eventid)
        readxh(ax=ax,filename=filename1,dirname=dirname,odirname=current_path,plot='True',align_index=9,phase='S',yy='intg',yy_index=6,ylabel='Distance (deg)',cut_b=50,cut_a=100,ofile='{}.pdf'.format(filename1),normalize='True',color_seism='blue')
        filename2 = "{}.BHT.dis.sflip.xh".format(eventid)
        readxh(ax=ax,filename=filename2,dirname=dirname,odirname=current_path,plot='True',align_index=3,phase='S',yy='intg',yy_index=6,ylabel='Distance (deg)',cut_b=50,cut_a=100,ofile='{}.pdf'.format(filename2),normalize='True',color_seism='black')
        fig.savefig('{}/{}.{}.pdf'.format(current_path,filename1,filename2.split('.')[3]))

#plot_overlap()
main()
#single_gridpoint()
