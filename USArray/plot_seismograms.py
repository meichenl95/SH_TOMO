#!/home/meichen/anaconda3/bin/python

def plot_seismograms_gcarc(**kwargs):

##-------------------------------##

# This function plot seismograms with y axis as distance(deg).

# Created by Meichen Liu on June 29th, 2019

##-------------------------------##

##parameters
# dirname		The directory where seismograms are stored
# filename		the name of files
# align			align to which phase/marked time.
#			'o' aligns to original time
#			'Tn' aligns to Tn time
# phase			phase to be aligned
# cut_b			seconds before the align time of the time window
#			(positive if before the align time)
# cut_a			seconds after the align time of the time window
# odirname		the directory to save the output file
# ofile			the name of output file
# normalize		True or False

    import numpy as np
    import matplotlib.pyplot as plt
    import obspy
    import os

    dirname = kwargs.get('dirname')
    filename = kwargs.get('filename')
    align = kwargs.get('align')
    phase = kwargs.get('phase')
    cut_b = kwargs.get('cut_b')
    cut_a = kwargs.get('cut_a')
    odirname = kwargs.get('odirname')
    ofile = kwargs.get('ofile')
    normalize = kwargs.get('normalize')
    
    os.chdir('{}'.format(dirname))

    fig,ax = plt.subplots(1,1,figsize=[5,8])
    st = obspy.read('{}'.format(filename))
    for tr in st:
        align_time = tr.stats.sac[align]
        start = tr.stats.sac['b']
        gcarc = tr.stats.sac['gcarc']
        delta = tr.stats.sac['delta']
        if align_time == '':
            print('No {} in the header file of {}'.format(align,tr.id))
            continue
        if start > align_time:
            print("Time window exceeds the range of {}.".format(tr.id))
            continue

        N1 = int((-1*start + align_time - cut_b)/delta)
        N2 = int((-1*start + align_time + cut_a)/delta)
        t = np.linspace(-1*cut_b,cut_a,num=len(tr.data[N1:N2]))
        if normalize == "True":
            if max(tr.data[N1:N2]) - min(tr.data[N1:N2]) != 0:
                ax.plot(t, tr.data[N1:N2]/(max(tr.data[N1:N2]) - min(tr.data[N1:N2]))*2 + gcarc, lw=0.4,c='k',alpha=0.5)
                ax.scatter(0,tr.data[int((-1*start+align_time)/delta)]/(max(tr.data[N1:N2]) - min(tr.data[N1:N2]))*2+gcarc,marker='o',color='red',s=0.5)
        elif normalize == "False":
            ax.plot(t,tr.data[N1:N2]*2*1e-4+gcarc,lw=0.4,c='k',alpha=0.3)
            ax.scatter(0,tr.data[int((-1*start+align_time)/delta)]*2*1e-4+gcarc,marker='o',color='red',s=0.5)
        else:
            print("Please specify normalize option: True or False!")

    ax.set_xlim([-1*cut_b-10,cut_a+10])
    ax.set_ylim([25,95])
    ax.set_xlabel('Time (s) aligned to {}'.format(phase))
    ax.set_ylabel('Distance (deg)')
    plt.savefig('{}/{}'.format(odirname,ofile))

def plot_seismograms_az(**kwargs):

##-------------------------------##

# This function plot seismograms with y axis as azimuth.

# Created by Meichen Liu on June 29th, 2019

##-------------------------------##

##parameters
# dirname		The directory where seismograms are stored
# filename		the name of files
# align			align to which phase/marked time.
#			'o' aligns to original time
#			'Tn' aligns to Tn time
# phase			phase to be aligned
# cut_b			seconds before the align time of the time window
#			(positive if before the align time)
# cut_a			seconds after the align time of the time window
# odirname		the directory to save the output file
# ofile			the name of output file
# normalize		True or False

    import numpy as np
    import matplotlib.pyplot as plt
    import obspy
    import os

    dirname = kwargs.get('dirname')
    filename = kwargs.get('filename')
    align = kwargs.get('align')
    phase = kwargs.get('phase')
    cut_b = kwargs.get('cut_b')
    cut_a = kwargs.get('cut_a')
    odirname = kwargs.get('odirname')
    ofile = kwargs.get('ofile')
    normalize = kwargs.get('normalize')
    
    os.chdir('{}'.format(dirname))

    fig,ax = plt.subplots(1,1,figsize=[5,8])
    st = obspy.read('{}'.format(filename))
    for tr in st:
        align_time = tr.stats.sac[align]
        start = tr.stats.sac['b']
        az = tr.stats.sac['az']
        delta = tr.stats.sac['delta']
        if align_time == '':
            print('No {} in the header file of {}'.format(align,tr.id))
            continue
        if start > align_time:
            print("Time window exceeds the range of {}.".format(tr.id))
            continue

        N1 = int((-1*start + align_time - cut_b)/delta)
        N2 = int((-1*start + align_time + cut_a)/delta)
        t = np.linspace(-1*cut_b,cut_a,num=len(tr.data[N1:N2]))
        if normalize == "True":
            if max(tr.data[N1:N2]) - min(tr.data[N1:N2]) != 0:
                ax.plot(t, tr.data[N1:N2]/(max(tr.data[N1:N2]) - min(tr.data[N1:N2]))*8 + az, lw=0.4,c='k',alpha=0.5)
                ax.scatter(0,tr.data[int((-1*start+align_time)/delta)]/(max(tr.data[N1:N2]) - min(tr.data[N1:N2]))*8+az,marker='o',color='red',s=0.5)
        elif normalize == "False":
            ax.plot(t,tr.data[N1:N2]*8*1e-4+az,lw=0.4,c='k',alpha=0.3)
            ax.scatter(0,tr.data[int((-1*start+align_time)/delta)]*8*1e-4+az,marker='o',color='red',s=0.5)
        else:
            print("Please specify normalize option: True or False!")

    ax.set_xlim([-1*cut_b-10,cut_a+10])
    ax.set_ylim([-10,370])
    ax.set_xlabel('Time (s) aligned to {}'.format(phase))
    ax.set_ylabel('Azimuth (deg)')
    plt.savefig('{}/{}'.format(odirname,ofile))
def main():

    import numpy as np
    import pandas as pd
    
    current_path = '/home/meichen/Research/SH_TOMO/USArray'
    path = '/home/meichen/work1/SH_TOMO/events'
    events_cat = pd.read_csv('{}/events_cat.txt'.format(current_path),skipinitialspace=True,header=None,sep=' ') 
    events_cat = np.array(events_cat)

    for eventid in events_cat[0:1,0]:
        plot_seismograms_gcarc(dirname='{}/event_{}/waveforms/SAC_files'.format(path,eventid),filename='*BHZ*.markP',align='t2',phase='P',cut_b=50,cut_a=200,odirname=current_path,ofile='event_{}_gcarc.pdf'.format(eventid),normalize='False')
        plot_seismograms_gcarc(dirname='{}/event_{}/waveforms/SAC_files'.format(path,eventid),filename='*BHZ*.markP',align='t2',phase='P',cut_b=50,cut_a=200,odirname=current_path,ofile='event_{}_gcarc_nor.pdf'.format(eventid),normalize='True')
        plot_seismograms_az(dirname='{}/event_{}/waveforms/SAC_files'.format(path,eventid),filename='*BHZ*.markP',align='t2',phase='P',cut_b=50,cut_a=200,odirname=current_path,ofile='event_{}_az.pdf'.format(eventid),normalize='False')
        plot_seismograms_az(dirname='{}/event_{}/waveforms/SAC_files'.format(path,eventid),filename='*BHZ*.markP',align='t2',phase='P',cut_b=50,cut_a=200,odirname=current_path,ofile='event_{}_az_nor.pdf'.format(eventid),normalize='True')

main()
