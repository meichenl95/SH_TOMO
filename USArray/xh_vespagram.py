#!/home/meichen/anaconda3/bin/python

def readvespaxh(**kwargs):

##---------------------------------------##

# This function read the xh file header.

# Created by Meichen Liu on July 6th, 2019.

##---------------------------------------##
##parameters
# filename	the name of the xh file
# dirname	the directory of filename
# odirname	the directory to save output file
# ofile		the name of output file

    import struct
    import matplotlib.pyplot as plt
    import numpy as np
    import os

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    odirname = kwargs.get('odirname')
    NPTS = kwargs.get('NPTS')
    samplingrate = kwargs.get('samplingrate')

    os.chdir('{}'.format(dirname))

    bihead_format = "f"*(NPTS+1+NPTS)
    t = np.arange(NPTS)*samplingrate -200

    with open("{}".format(filename),"rb") as fp:
        i = 0
        while True:
            bihead = fp.read(4*NPTS+4+4*NPTS)
            if bihead:
                head = struct.unpack(bihead_format,bihead)
                slowness = np.float(head[0])
                env = np.log10(head[1:NPTS+1])
                stack = (-1)*np.sign(head[1+NPTS::])*np.log10(np.abs(head[NPTS+1::]))

                x = np.vstack((t,np.ones(len(t))*slowness,env)).T
                y = np.vstack((t,np.ones(len(t))*slowness,stack)).T
                if(i == 0):
                    b = x
                    c = y
                else:
                    b = np.concatenate((b,x),axis=0)
                    c = np.concatenate((c,y),axis=0)
                i = i +1
            else:
                break
    np.savetxt('gmtvespa_env.txt',b)
    np.savetxt('gmtvespa_noenv.txt',c)
    fp.close()

def main():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    path = '/home/meichen/work1/SH_TOMO/Combined_xh/vesp_rp_all'
    current_path = '/home/meichen/Research/SH_TOMO/USArray'
    
    dirname = '{}'.format(path)
    filename = "vesp.anchor_85.dis"
    readvespaxh(filename=filename,dirname=dirname,odirname=current_path,NPTS=7000,samplingrate=0.1)

main()
