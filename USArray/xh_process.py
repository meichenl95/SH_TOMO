#!/home/meichen/anaconda3/bin/python

def mk_path_rmexist(path):

##---------------------------------##

# This function check if path exist. If so, delete it and make a new directory.

##---------------------------------##

    import os
    import subprocess
    isExist = os.path.exists(path)
    if isExist:
        subprocess.call(['rm -r {}'.format(path)],shell=True)
    os.makedirs(path)

def mk_path_retainexist(path):

##---------------------------------##

# This function check if path exist. If not, create the directory
##---------------------------------##
    import os
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)

def xh_process_stf(**kwargs):

##---------------------------------##

# This function write sac_2xh commands into a .sh file with corresponding names.

# Created by Meichen Liu on June 18th, 2019

##---------------------------------##

##parameters
# cp_from_path		the directory to copy xh files from
# dirname 		the directory to save xh files
# seis_type		the type of seismograms
# eventid		the id of the event
# max_takeoff_TP	the maximum accepted angle between takoff and min(T,P)
# min_takeoff_nodal	the minimum accepted angle between takeoff and nodal planes
# sum_ratio		the maximum ratio between the sum of noise and signal
# max_ratio		the maximum ratio between the max values of noise and signal
# num_seism		the minimum number of seismograms when stacking for source-time-function
# events_selected		the file to save events satisfying selection criteria

    import os
    import subprocess
    import glob

    cp_from_path = kwargs.get('cp_from_path')
    dirname = kwargs.get('dirname')
    eventid = kwargs.get('eventid')
    chn = kwargs.get('chn')
    seis_type = kwargs.get('seis_type')
    max_takeoff_TP = kwargs.get('max_takeoff_TP')
    min_takeoff_nodal = kwargs.get('min_takeoff_nodal')
    sum_ratio = kwargs.get('sum_ratio')
    max_ratio = kwargs.get('max_ratio')
    num_seism = kwargs.get('num_seism')
    events_selected = kwargs.get('events_selected')

    mk_path_rmexist('{}'.format(dirname))
    os.chdir('{}'.format(dirname))
        
    f = open('xh_process.sh','w')
    f.write("#!/bin/bash\n")
   
    xh_id = '{}.{}.{}'.format(eventid,chn,seis_type)

    # copy original xh files to new directory
    f.write("cp {}/{}.xh .\n".format(cp_from_path,xh_id))

    # save the angle between takeoff and TBP axis to flt[0], flt[1], and flt[2]
    # save the angle between takeoff and nodal planes to flt[3] and flt[4]
    # save takeoff angle to flt[5]
    # xh_takeoff_tbp save distance and azimuth to flt[6] and flt[7]
    f.write("xh_takeoff_tbp {}.xh P {}.takeoff_tbp.xh {}\n".format(xh_id,xh_id,eventid))
    
    # pick the peak time around P(tpcks[0]) and S(tpcks[1]) wave by envelope and save to tpcks[2] and tpcks[3]
    f.write("xh_envelope_pick {}.takeoff_tbp.xh {}.peakP.xh 0 2 10 30\n".format(xh_id,xh_id))
    f.write("xh_envelope_pick {}.peakP.xh {}.peakS.xh 1 3 10 30\n".format(xh_id,xh_id))

    # flip traces according to tpcks[N]
    f.write("xh_flip {}.peakS.xh {}.pflip.xh 2\n".format(xh_id,xh_id))

    # Signal-to-noise ratio to tpcks[2]
    f.write("xh_snr {}.pflip.xh {}.snrP.xh 2 5 5 -60 30\n".format(xh_id,xh_id))

    # xh_select_P select traces satisfy criteria.
    # We did signal-noise ratio above, so we do not snr anymore here.
    f.write("xh_select_P {}.snrP.xh {}.select_P.xh {} {} {} {}\n".format(xh_id,xh_id,max_takeoff_TP,min_takeoff_nodal,sum_ratio,max_ratio))

    # xh_boxcar extracted P wave
    f.write("xh_boxcar {}.select_P.xh {}.boxcar.xh 2 10 10 \n".format(xh_id,xh_id))

    # xh_normalize normalize extracted P waves
    f.write("xh_normalize {}.boxcar.xh {}.nor.xh \n".format(xh_id,xh_id))

    # xh_header gives the number of seismograms stacked
    # we require at least 10 seismograms
    # save the number of seismograms to num_seism
    # if begins
    f.write("num_seism=`xh_header {}.nor.xh` \n".format(xh_id))
    f.write("if [ $num_seism -ge {} ];then\n".format(num_seism))
    f.write("echo {} $num_seism >> {}\n".format(eventid,events_selected))

    # xh_stack average aligned, normalized P waves to get stf
    f.write("xh_stack {}.nor.xh {}.stack 0 2 15 15 \n".format(xh_id,xh_id))

    # xh_deconv deconvolve source-time-function from seismograms
    f.write("xh_deconv {}.select_P.xh {}.stack {}.deconv.xh 0.0001 \n".format(xh_id,xh_id,xh_id))

    f.close()
        
def xh_process_deconv(**kwargs):

##---------------------------------##

# This function write sac_2xh commands into a .sh file with corresponding names.

# Created by Meichen Liu on June 18th, 2019

##---------------------------------##

##parameters
# cp_from_path		the directory to copy xh files from
# dirname 		the directory to save xh files
# seis_type		the type of seismograms
# eventid		the id of the event
# events_selected		the file to save events satisfying selection criteria

    import os
    import subprocess
    import glob

    cp_from_path = kwargs.get('cp_from_path')
    dirname = kwargs.get('dirname')
    eventid = kwargs.get('eventid')
    chn = kwargs.get('chn')
    seis_type = kwargs.get('seis_type')
    events_selected = kwargs.get('events_selected')

    os.chdir('{}'.format(dirname))
        
    f = open('xh_process.sh','w')
    f.write('#!/bin/bash\n')
   
    xh_id = '{}.{}.{}'.format(eventid,chn,seis_type)

    # copy original xh files to new directory
    f.write("rm *.xh\n")
    f.write("cp {}/{}.xh .\n".format(cp_from_path,xh_id))
    
    # mark travel time of S to tpcks[1]
    # P to tpcks[0]
    # S^660S to tpcks[5]
    # S^410S to tpcks[6]
    # sS to tpcks[7]
    # ScS to tpcks[8]
    # sScS to tpcks[4]
    # SS to tpcks[9]
    # SSv410s to tpcks[10]
    # SSv660s to tpcks[11]
#    f.write("xh_marktt {}.xh {}.markP.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.P 0 0 5\n".format(xh_id,xh_id))
    f.write("xh_marktt {}.xh {}.markS.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.S 1 0 5\n".format(xh_id,xh_id))
#    f.write("xh_marktt {}.markS.xh {}.markS660S.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.S^660S 5 1 10\n".format(xh_id,xh_id))
#    f.write("xh_marktt {}.markS660S.xh {}.markS410S.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.S^410S 6 1 10\n".format(xh_id,xh_id))
#    f.write("xh_marktt {}.markS410S.xh {}.marksS.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.sS 7 1 10\n".format(xh_id,xh_id))
#    f.write("xh_marktt {}.marksS.xh {}.markSS.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.SS 9 1 10\n".format(xh_id,xh_id))
#    f.write("xh_marktt {}.markSS.xh {}.markScS.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.ScS 8 1 10\n".format(xh_id,xh_id))
#    f.write("xh_marktt {}.markScS.xh {}.marksScS.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.sScS 4 1 10\n".format(xh_id,xh_id))
#    f.write("xh_marktt {}.marksScS.xh {}.markSSv410s.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.SSv410s 10 1 5\n".format(xh_id,xh_id))
#    f.write("xh_marktt {}.markSSv410s.xh {}.markSSv660s.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.SSv660s 11 1 5\n".format(xh_id,xh_id))
#    f.write("rm {}.markP.xh {}.markS.xh {}.markS660S.xh {}.markS410S.xh {}.marksS.xh {}.markSS.xh {}.markScS.xh {}.marksScS.xh {}.markSSv410s.xh\n".format(xh_id,xh_id,xh_id,xh_id,xh_id,xh_id,xh_id,xh_id,xh_id))
    # save the angle between takeoff and TBP axis to flt[0], flt[1], and flt[2]
    # save the angle between takeoff and nodal planes to flt[3] and flt[4]
    # save takeoff angle to flt[5]
    # xh_takeoff_tbp save distance and azimuth to flt[6] and flt[7]
    # xh_takeoff_tbp only save station distance>50 degrees
    # xh_takeoff_tbp select 25.0<slat<50.0 -130.0<slon<-65.0
    # xh_takeoff_tbp also save "my evtcode" to XH file.
    f.write("xh_takeoff_tbp {}.markS.xh P {}.takeoff_tbp.xh {}\n".format(xh_id,xh_id,eventid))
#    f.write("rm {}.markS.xh\n".format(xh_id))

    # xh_changewavf change h.wavf 
    f.write("xh_changewavf {}.takeoff_tbp.xh {}.changewavf.xh dis\n".format(xh_id,xh_id))
    f.write("rm {}.takeoff_tbp.xh\n".format(xh_id))

    # xh_bpfilter bandpass to desired freq range
    # According to experiment on synthetic waveforms, there is no difference in resulting waveforms in terms of the order of xh_bpfilter and xh_disp2vel
    f.write("xh_bpfilter {}.changewavf.xh {}.bp.xh -t 10 -f 10 20 80 120\n".format(xh_id,xh_id))

    # pick the absolute peak time around displacement S(tpcks[1]) and save to tpcks[3]
    f.write("xh_peak {}.bp.xh {}.peakS.xh 1 3 10 40 0 1\n".format(xh_id,xh_id))
#    f.write("xh_peak {}.bp.xh {}.peakS.xh 1 3 50 50 0 0.7\n".format(xh_id,xh_id))
    f.write("rm {}.bp.xh\n".format(xh_id))

##    # normalize traces according to tpcks[N]
##    f.write("xh_normalize {}.peakS.xh {}.norm.xh 3\n".format(xh_id,xh_id))
##    f.write("rm {}.peakS.xh\n".format(xh_id))
##
##    # displacement to velocity
##    f.write("xh_disp2vel {}.norm.xh {}.{}.vel.norm.xh -t 10\n".format(xh_id,eventid,chn))
##    f.write("rm {}.norm.xh\n".format(xh_id))
##
##    # pick the peak time in velocity
##    f.write("xh_peak_vel {}.{}.vel.norm.xh {}.{}.vel.peakS.xh 3 3 30 0 1\n".format(eventid,chn,eventid,chn))
##    f.write("rm {}.{}.vel.norm.xh\n".format(eventid,chn))

    # normalize before doing root mean squre
    f.write("xh_normalize {}.peakS.xh {}.norm.xh 3\n".format(xh_id,xh_id))

    # xh_snr select high snr traces
    f.write("xh_snr_rms {}.norm.xh {}.snrS.xh 3 5 5 -150 100 5 2.5 50 100\n".format(xh_id,xh_id))

    f.close()
    subprocess.call(['sh xh_process.sh'],shell=True)
#    subprocess.call(['rm xh_process.sh'],shell=True)
        
def combine_xh(**kwargs):
## parameters
# path		the directory of the file to read in
# filename	the name of the file to read in
# dirname	the directory where to create and save .sh file
# xhappend	the appendix of the xh files that are involved
# save_path	the directory of the xh files where they are saved
# num_seism	the minimum number of seismograms for one event

    import os
    import pandas as pd
    import numpy as np
    import subprocess

    path = kwargs.get('path')
    save_path = kwargs.get('save_path')
    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    xhappend = kwargs.get('xhappend')
    num_seism = kwargs.get('num_seism')

    mk_path_retainexist('{}'.format(dirname))
    os.chdir('{}'.format(dirname))

    events_selected = pd.read_csv('{}/{}.txt'.format(path,filename),skipinitialspace=True,header=None,sep=' ')
    events_selected = np.array(events_selected)
    # One event should have at least num_seism+1 records
    events_selected = events_selected[events_selected[:,1]>num_seism]

    # Write a bash script
    f = open('add.sh','w')
    f.write('#!/bin/bash\n')
    f.write('ulimit -S -s unlimited\n')
    f.write('ulimit -S -n 1048576\n')
    for eventid in events_selected[:,0]:
        f.write('cp {}/event_{}/waveforms/XH_files/{}.{} {}.{}.temp \n'.format(save_path,eventid,eventid,xhappend,eventid,xhappend))

    f.write('cat *.{}.temp > all.{}\n'.format(xhappend,xhappend))
    f.write('rm *.{}.temp\n'.format(xhappend))
    f.write('\n')
    f.close()
    subprocess.call(['sh add.sh'],shell=True)
#    subprocess.call(['rm add.sh'],shell=True)
    

def main():

    import sys
    import pandas as pd
    import numpy as np
    import subprocess

    path = '/home/meichen/Research/SH_TOMO/USArray'
    save_path = '/home/meichen/work1/SH_TOMO/events'
    events_selected = '{}/events_selected.txt'.format(path)

#    # save selected events
#    subprocess.call(['echo -n > {}'.format(events_selected)],shell=True)

#    events_cat = pd.read_csv('{}/events_cat/events_cat_300_800_2010_2018_Mw6075.txt'.format(path),skipinitialspace=True,header=None,sep=' ')
#    events_cat = np.array(events_cat)
#    for eventid in events_cat[:,0]:
#        folder_path = '{}/event_{}/waveforms'.format(save_path,eventid)
#        seis_type = 'dis'
#        xh_process_deconv(cp_from_path='{}/SAC_files'.format(folder_path),dirname='{}/XH_files'.format(folder_path),eventid=eventid,chn='BHT',seis_type=seis_type,events_selected=events_selected)

    # Adding XH files
    combine_xh(path=path,save_path=save_path,filename='events_selected',dirname='/home/meichen/work1/SH_TOMO/Combined_xh',xhappend='BHT.dis.snrS.xh',num_seism=20)

main()
