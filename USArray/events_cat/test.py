#!/home/meichen/anaconda3/bin/python

import numpy as np
import pandas as pd

f = open('run.sh','w')
f.write('#!/bin/bash\n')

data = pd.read_csv('events_cat_all_5560.txt',skipinitialspace=True,sep=' ',header=None)
f.write('cd /home/meichen/work1/SH_TOMO/events\n')

for filename in data[0]:
    f.write('rm event_{}\n'.format(filename))

f.close()
