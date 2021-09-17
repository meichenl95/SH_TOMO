#!/home/meichen/anaconda3/bin/python

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

    dirname = kwargs.get('dirname')
    modelid = kwargs.get('modelid')
    timeS = kwargs.get('timeS')
    component = kwargs.get('component')

    os.chdir('{}'.format(dirname))
   
    f = open('xh_process.sh','w')
    f.write('#!/bin/bash\n')

    # rotate ZNE --> ZRT
    f.write("xh_rotate {}.dis.xh {}.rotate.dis.xh\n".format(modelid,modelid))
    # grep T component
    f.write("xh_grepchan {}.rotate.dis.xh {}.{}.dis.xh -c {}\n".format(modelid,modelid,component,component))
    
    # mark travel time of S to tpcks[1]
    # S^660S to tpcks[5]
    # S^410S to tpcks[6]
    # sS to tpcks[7]
    # ScS to tpcks[8]
    # SS to tpcks[9]
    # SSv410s to tpcks[10]
    # SSv660s to tpcks[11]
    f.write("xh_tpck {}.{}.dis.xh ttc_file_S40RTS_S {}.markS.{}.dis.xh {}\n".format(modelid,component,modelid,component,timeS))
    f.write("xh_marktt {}.markS.{}.dis.xh {}.markS660S.{}.dis.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.S^660S 5 1 10\n".format(modelid,component,modelid,component))
    f.write("xh_marktt {}.markS660S.{}.dis.xh {}.markS410S.{}.dis.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.S^410S 6 1 10\n".format(modelid,component,modelid,component))
    f.write("xh_marktt {}.markS410S.{}.dis.xh {}.marksS.{}.dis.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.sS 7 1 10\n".format(modelid,component,modelid,component))
    f.write("xh_marktt {}.marksS.{}.dis.xh {}.markSS.{}.dis.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.SS 9 1 10\n".format(modelid,component,modelid,component))
    f.write("xh_marktt {}.markSS.{}.dis.xh {}.markScS.{}.dis.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.ScS 8 1 10\n".format(modelid,component,modelid,component))
    f.write("xh_marktt {}.markScS.{}.dis.xh {}.marksScS.{}.dis.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem.sScS 4 1 10\n".format(modelid,component,modelid,component))
    f.write("xh_marktt {}.marksScS.{}.dis.xh {}.markSSv410s.{}.dis.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.SSv410s 10 1 5\n".format(modelid,component,modelid,component))
    f.write("xh_marktt {}.markSSv410s.{}.dis.xh {}.markSSv660s.{}.dis.xh ~/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.SSv660s 11 1 5\n".format(modelid,component,modelid,component))

    # xh_bpfilter bandpass to desired freq range
    f.write("xh_bpfilter {}.markSSv660s.{}.dis.xh {}.bp.{}.dis.xh -t 10 -f 10 20 80 120\n".format(modelid,component,modelid,component))

    # pick the absolute peak time around S(tpcks[1]) wave and save to tpcks[3]
    f.write("xh_peak {}.bp.{}.dis.xh {}.peakS.{}.dis.xh 1 3 10 30 0 1\n".format(modelid,component,modelid,component))

    # normalize traces according to tpcks[N]
    f.write("xh_normalize {}.peakS.{}.dis.xh {}.norm.{}.dis.xh 3\n".format(modelid,component,modelid,component))

    # displacement to velocity
    f.write("xh_disp2vel {}.norm.{}.dis.xh {}.norm.{}.vel.xh -t 10\n".format(modelid,component,modelid,component))

    # pick the peak time in velocity
    f.write("xh_peak_vel {}.norm.{}.vel.xh {}.peakS.{}.vel.xh 3 3 20 0 1\n".format(modelid,component,modelid,component))

    # normalize to tpck[3] before doing snr
    f.write("xh_normalize {}.peakS.{}.vel.xh {}.norm.{}.vel.xh 3\n".format(modelid,component,modelid,component))

    # xh_snr select high snr traces
    f.write("xh_snr_rms {}.norm.{}.vel.xh {}.snrS.{}.vel.xh 3 6 6 -130 100 6 3 30 200\n".format(modelid,component,modelid,component))

    f.close()
#    subprocess.call(['sh xh_process.sh'],shell=True)
        
def main():

    import sys
    import pandas as pd
    import numpy as np
    import subprocess

    model='MTZ_8X8_NOCRUST_10s'
    modelfln='mtz_8x8_nocrust_10s'
    timeS='theo'
    component='BXT'
    save_path = '/home/meichen/work1/SH_TOMO/synthetics/XH_files'

    xh_process_deconv(dirname='{}/{}'.format(save_path,model),modelid=modelfln,timeS=timeS,component=component)

main()
