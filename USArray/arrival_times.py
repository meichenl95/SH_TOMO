#!/home/meichen/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt

path = "/home/meichen/Utils/TauP_calculations/TravelTimeTables"
file_S = "tt.prem.S"
file_ScS = "tt.prem.ScS"
file_sScS = "tt.prem.sScS"
file_SS = "tt.prem.SS"
file_sS = "tt.prem.sS"
file_rever = "rt.prem_ml..s"
file_S410S = "tt.prem_ml.S^410S"
file_S660S = "tt.prem_ml.S^660S"
file_S150S = "tt.prem_ml.S^150S"
file_SSv410sSv410s = "tt.prem_ml.SSv410sSv410s"
file_PSSv660s = "tt.prem_ml.PSSv660s"
edep = 20

# start plot
fig,ax = plt.subplots(1,1,figsize=[7,11])

# SSv410, SSv660
depth1 = 410
depth2 = 660
depth3 = 150
rt = np.genfromtxt('{}/{}'.format(path,file_rever))
reflection_depth = rt[:,0]
event_depth = rt[:,1]
distance = rt[:,2]
time_resp_S = rt[:,4]
event_depth_410 = event_depth[np.abs(reflection_depth-depth1)<0.1]      
time_resp_S_410 = time_resp_S[np.abs(reflection_depth-depth1)<0.1]
distance_410 = distance[np.abs(reflection_depth-depth1)<0.1]
time_resp_S_410 = time_resp_S_410[np.abs(event_depth_410-edep)<0.1]
distance_410 = distance_410[np.abs(event_depth_410-edep)<0.1]
distance_410 = distance_410[time_resp_S_410>0.0]
time_resp_S_410 = time_resp_S_410[time_resp_S_410>0.0]
ax.plot(time_resp_S_410,distance_410,color='k',lw=1)
event_depth_660 = event_depth[np.abs(reflection_depth-depth2)<0.1]
time_resp_S_660 = time_resp_S[np.abs(reflection_depth-depth2)<0.1]
distance_660 = distance[np.abs(reflection_depth-depth2)<0.1]
time_resp_S_660 = time_resp_S_660[np.abs(event_depth_660-edep)<0.1]
distance_660 = distance_660[np.abs(event_depth_660-edep)<0.1]
distance_660 = distance_660[time_resp_S_660>0.0]
time_resp_S_660 = time_resp_S_660[time_resp_S_660>0.0]
ax.plot(time_resp_S_660,distance_660,color='k',lw=1)
event_depth_150 = event_depth[np.abs(reflection_depth-depth3)<0.1]
time_resp_S_150 = time_resp_S[np.abs(reflection_depth-depth3)<0.1]
distance_150 = distance[np.abs(reflection_depth-depth3)<0.1]
time_resp_S_150 = time_resp_S_150[np.abs(event_depth_150-edep)<0.1]
distance_150 = distance_150[np.abs(event_depth_150-edep)<0.1]
distance_150 = distance_150[time_resp_S_150>0.0]
time_resp_S_150 = time_resp_S_150[time_resp_S_150>0.0]
ax.plot(time_resp_S_150,distance_150,color='k',lw=1)

# S
S = np.genfromtxt('{}/{}'.format(path,file_S))
tt_S = S[:,2]
depth_S = S[:,1]
tt_S = tt_S[np.abs(depth_S-edep)<0.1]
ax.plot([0,0],[0,180],'r-')

# SS
SS = np.genfromtxt('{}/{}'.format(path,file_SS))
tt = SS[:,2]
depth = SS[:,1]
distance = SS[:,0]
tt = tt[np.abs(depth-edep)<0.1]
distance = distance[np.abs(depth-edep)<0.1]
tt_S_temp = tt_S[tt<9999]
distance = distance[tt<9999]
tt = tt[tt<9999]
tt = tt[tt_S_temp>0]
distance = distance[tt_S_temp>0]
tt_S_temp = tt_S_temp[tt_S_temp>0]
ax.plot(tt-tt_S_temp,distance,color='k',lw=1)

# SSv410sSv410s
SSv410sSv410s = np.genfromtxt('{}/{}'.format(path,file_SSv410sSv410s))
tt = SSv410sSv410s[:,2]
depth = SSv410sSv410s[:,1]
distance = SSv410sSv410s[:,0]
tt = tt[np.abs(depth-edep)<0.1]
distance = distance[np.abs(depth-edep)<0.1]
tt_S_temp = tt_S[tt<9999]
distance = distance[tt<9999]
tt = tt[tt<9999]
tt = tt[tt_S_temp>0]
distance = distance[tt_S_temp>0]
tt_S_temp = tt_S_temp[tt_S_temp>0]
ax.plot(tt-tt_S_temp,distance,color='k',lw=1)

# sS
sS = np.genfromtxt('{}/{}'.format(path,file_sS))
tt = sS[:,2]
depth = sS[:,1]
distance = sS[:,0]
tt = tt[np.abs(depth-edep)<0.1]
distance = distance[np.abs(depth-edep)<0.1]
tt_S_temp = tt_S[tt<9999]
distance = distance[tt<9999]
tt = tt[tt<9999]
tt = tt[tt_S_temp>0]
distance = distance[tt_S_temp>0]
tt_S_temp = tt_S_temp[tt_S_temp>0]
ax.plot(tt-tt_S_temp,distance,color='k',lw=1)

# sScS
sScS = np.genfromtxt('{}/{}'.format(path,file_sScS))
tt = sScS[:,2]
depth = sScS[:,1]
distance = sScS[:,0]
tt = tt[np.abs(depth-edep)<0.1]
distance = distance[np.abs(depth-edep)<0.1]
tt_S_temp = tt_S[tt<9999]
distance = distance[tt<9999]
tt = tt[tt<9999]
tt = tt[tt_S_temp>0]
distance = distance[tt_S_temp>0]
tt_S_temp = tt_S_temp[tt_S_temp>0]
ax.plot(tt-tt_S_temp,distance,color='k',lw=1)

# S^660S
S660S = np.genfromtxt('{}/{}'.format(path,file_S660S))
tt = S660S[:,2]
depth = S660S[:,1]
distance = S660S[:,0]
tt = tt[np.abs(depth-edep)<0.1]
distance = distance[np.abs(depth-edep)<0.1]
tt_S_temp = tt_S[tt<9999]
distance = distance[tt<9999]
tt = tt[tt<9999]
tt = tt[tt_S_temp>0]
distance = distance[tt_S_temp>0]
tt_S_temp = tt_S_temp[tt_S_temp>0]
ax.plot(tt-tt_S_temp,distance,color='k',lw=1)

# S^410S
S410S = np.genfromtxt('{}/{}'.format(path,file_S410S))
tt = S410S[:,2]
depth = S410S[:,1]
distance = S410S[:,0]
tt = tt[np.abs(depth-edep)<0.1]
distance = distance[np.abs(depth-edep)<0.1]
tt_S_temp = tt_S[tt<9999]
distance = distance[tt<9999]
tt = tt[tt<9999]
tt = tt[tt_S_temp>0]
distance = distance[tt_S_temp>0]
tt_S_temp = tt_S_temp[tt_S_temp>0]
ax.plot(tt-tt_S_temp,distance,color='k',lw=1)

# S^150S
S150S = np.genfromtxt('{}/{}'.format(path,file_S150S))
tt = S150S[:,2]
depth = S150S[:,1]
distance = S150S[:,0]
tt = tt[np.abs(depth-edep)<0.1]
distance = distance[np.abs(depth-edep)<0.1]
tt_S_temp = tt_S[tt<9999]
distance = distance[tt<9999]
tt = tt[tt<9999]
tt = tt[tt_S_temp>0]
distance = distance[tt_S_temp>0]
tt_S_temp = tt_S_temp[tt_S_temp>0]
ax.plot(tt-tt_S_temp,distance,color='k',lw=1)

ax.set_xlim([-50,400])
ax.set_ylim([110,60])
ax.set_xlabel('Time (s)',size=14)
ax.set_ylabel('Distance',size=14)
ax.text(-20,108,'S',fontsize=8)
ax.text(20,100,'sS',fontsize=8)
ax.text(28,88,'sScS',fontsize=8)
ax.text(22,70,'ScS',fontsize=8)
ax.text(105,90,'S660S',fontsize=8)
ax.text(340,105,'S410S',fontsize=8)
ax.text(290,95,'S150S',fontsize=8)
ax.text(130,103,'SSv410s',fontsize=8)
ax.text(260,108,'SSv660s',fontsize=8)
ax.text(70,100,'SSv150s',fontsize=8)
ax.text(320,80,'SS',fontsize=8)
ax.tick_params('x',bottom=True,top=True)
ax.tick_params('y',bottom=True,top=True)
plt.savefig('arrival_times.pdf')
