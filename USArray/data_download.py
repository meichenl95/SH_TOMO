#!/home/meichen/anaconda3/bin/python

def mk_path(path):

##**************************************##

# This function delete the original directory and make a new one.

##**************************************##
    import os
    import subprocess
    isExist = os.path.exists(path)
    if isExist:
        subprocess.call(['rm -r {}'.format(path)],shell=True)
    os.makedirs(path)

def mk_path_exist(path):

##**************************************##

# This function make a new directory if it does not exist.

##**************************************##
    import os
    import subprocess
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)

def download_data(**kwargs):

##**************************************##
# Download data through fdsn
# Original download data is miniseed and stationxml files, then change to 
# standard SAC files. Channels with gaps and coverage smaller than 95% of desired
# time range are excluded. Minimum interstation distance is 100 meters.

# Created by Meichen Liu on June 8th, 2019

##**************************************##

## parameters
# event_info	in array, [event time, event latitude, event longitude]
# shape		shape of station domain
#		rec	rectangular domain
#		cir	circular domain
# minlat	minimum latitude if shape=rec
# maxlat	maximum latitude if shape=rec
# minlon	minimum longitude if shape=rec
# maxlon	maximum longitude if shape=rec
# cenlat	central latitude if shape=cir
# cenlon	central longitude if shape=cir
# minrad	minimum radius in degrees if shape=cir
# maxrad	maximum radius in degrees if shape=cir
# start		seconds before event time
# end		seconds after event time
# channel	the channel choose for download,"BH[NEZ]","HH[NEZ]",see IRIS 
#		classification
# network	the network code
# client	All currently availabel web service providers are:
#		"BGR", "EMSC", "ETH", "GEONET", "GFZ", "ICGC", "INGV", "IPGP",
#		"IRIS", "ISC", "KNMI", "KOERI", "LMU", "NCEDC", "NIEP", "NOA",
#		"ODC", "ORFEUS", "RESIF", "SCEDC", "TEXNET", "USGS", "USP"
# save_path	the path to save files including miniseed, stationxml, and SAC files.

    import obspy
    import numpy as np
    from obspy.clients.fdsn.mass_downloader import CircularDomain, RectangularDomain, Restrictions, MassDownloader
    from sys import argv

    event_info = kwargs.get('event_info')
    shape = kwargs.get('shape')
    if shape == 'rec':
        minlat = kwargs.get('minlat')
        maxlat = kwargs.get('maxlat')
        minlon = kwargs.get('minlon')
        maxlon = kwargs.get('maxlon')
    elif shape == 'cir':
        cenlat = kwargs.get('cenlat')
        cenlon = kwargs.get('cenlon')
        minrad = kwargs.get('minrad')
        maxrad = kwargs.get('maxrad')
    start = kwargs.get('start')
    end = kwargs.get('end')
    channel = kwargs.get('channel')
    client = kwargs.get('client')
    save_path = kwargs.get('save_path')
    mk_path(save_path)
#    network = ['5A','6E','7A','7C','AE','AG','AK','AR','AT','AV','AZ','BK','CI','CN','EM','ET','II','IM','IU','IW','LB','LD','LM','MU','N4','NN','NY','OH','OK','PO','SC','TA','TX','UO','US','UU','UW','WY','X8','XA','XD','XE','XG','XI','XN','XO','XQ','XR','XT','XU','XV','YE','YG','YH','YO','YQ','YT','YW','YX','Z9','ZE','ZG','ZH','ZI','ZK','ZL','ZZ']
    network = ['_US-TA']
        
    origin_time = obspy.UTCDateTime(event_info[0])
    event_lat = event_info[1]
    event_lon = event_info[2]
    if shape == 'rec':
        domain = RectangularDomain(minlatitude=minlat,maxlatitude=maxlat,minlongitude=minlon,maxlongitude=maxlon)
    elif shape == 'cir':
        domain = CircularDomain(latitude=cenlat,longitude=cenlon,minradius=minrad,maxradius=maxrad)
   
    for i in np.arange(len(network)):
        net = network[i]
        restrictions = Restrictions(starttime=origin_time-start,endtime=origin_time+end,reject_channels_with_gaps=True,minimum_length=0.95,minimum_interstation_distance_in_m=1E2,network=net,channel_priorities=["{}".format(channel)],location_priorities=["","00","10","20","01","02"])
        mdl=MassDownloader(providers=['{}'.format(client)])
        mdl.download(domain,restrictions,mseed_storage="{}/waveforms".format(save_path),stationxml_storage="{}/stations".format(save_path))

def mseed2sac(**kwargs):

##**************************************##

# This function transform miniseed files to SAC files.
# First the javascript "stationxml-seed-converter-2.0.0.jar" produces dataless 
# files from stationxml files. Then "rdseed" using dataless files to tranform 
# mseed to SAC. Make sure rdseed is installed. Note that SAC files created from
# mseed lack some event information including event latitude, event longitude, 
# and event depth.

# Created by Meichen Liu on June 8th, 2019

##**************************************##

##parameters
# save_path	same as save_path in the function download_data
# jar_path	where the script "stationxml-seed-converter-2.0.0.jar is stored

    import subprocess
    import glob
    import os

    save_path = kwargs.get('save_path')
    jar_path = kwargs.get('jar_path')
    os.chdir('{}'.format(save_path))
    mk_path('{}/waveforms/SAC_files'.format(save_path))
    # change directory to save_path because rdseed does not work with home directory in the path

    for stnxml in glob.glob('{}/stations/*'.format(save_path)):
        stationname = stnxml.split('/')[-1]
        nw = stationname.split('.')[0]
        stn = stationname.split('.')[1]
        subprocess.call(['java','-jar','{}/stationxml-seed-converter-2.0.4-SNAPSHOT.jar'.format(jar_path),'--input','{}/stations/{}.{}.xml'.format(save_path,nw,stn),'--output','{}/waveforms/{}.{}.dataless'.format(save_path,nw,stn)])


        for filename in glob.glob('{}/waveforms/{}.{}.*.mseed'.format(save_path,nw,stn)):
            mseedfile = filename.split('/')[-1]
            subprocess.call(['rdseed','-pRdf','waveforms/{}'.format(mseedfile),'-z','1','-g','waveforms/{}.{}.dataless'.format(nw,stn),'-q','waveforms/SAC_files'])

def rename(**kwargs):

##**************************************##

# This function rename sac files to simpler, clearer filenames.

# Created by Meichen Liu on June 12nd, 2019 based on sac-manualv3.6

##**************************************##

##parameters
# dirname	the directory SAC files are saved

    import os
    import sys
    import glob

    dirname = kwargs.get('dirname')

    os.putenv("SAC_DISPLAY_COPYRIGHT","O")
    os.chdir('{}'.format(dirname))

    for filename in glob.glob("*.SAC"):
        nw, stn, loc, chn = filename.split('.')[6:10]
        os.rename(filename,"%s.%s.%s.%s.SAC.raw" % (nw, stn, loc, chn))


def add_info(**kwargs):

##**************************************##

# This function add event info to SAC files, including event latitude, event
# longitude, event depth, event magnitude, event original time.

# Created by Meichen Liu on June 10th, 2019 based on sac-manual-v3.6
##**************************************##

##parameters
# filename	the name of the SAC file
# dirname	the directory where SAC files are saved
# origin	the original time of the event
# evlo		the longitude of the event
# evla		the latitude of the event
# mag		the magnitude of the event

    import os
    import sys
    import datetime
    import subprocess
    import glob
    import numpy as np

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    origin = kwargs.get('origin')
    evlo = kwargs.get('evlo')
    evla = kwargs.get('evla')
    evdp = kwargs.get('evdp')
    mag = kwargs.get('mag')

    os.putenv("SAC_DISPLAY_COPYRIGHT","O")
    os.chdir('{}'.format(dirname))

    o = datetime.datetime.strptime(origin,'%Y-%m-%dT%H:%M:%S')
    # calculate which day in a year is the occurence date
    jday = o.strftime("%j")

    p = subprocess.Popen(['sac'],stdin=subprocess.PIPE)
    s = "wild echo off \n"

    filelist = glob.glob('{}'.format(filename))
    if len(filelist)>20:
        for i in np.arange(round(len(filelist)/20)):
            s += "r %s \n" % filelist[i*20:i*20+20]
            s += "synchronize \n"
            s += "ch o gmt %s %s %s %s %s \n" % (o.year, jday, o.hour, o.minute, o.second)
            s += "ch allt (0 - &1,o&) iztype IO \n"
            s += "ch evlo %s evla %s evdp %s mag %s \n" % (evlo, evla, evdp, mag)
            s += "wh \n"
        s += "r %s \n" % filelist[i*20+20::]
        s += "synchronize \n"
        s += "ch o gmt %s %s %s %s %s \n" % (o.year, jday, o.hour, o.minute, o.second)
        s += "ch allt (0 - &1,o&) iztype IO \n"
        s += "ch evlo %s evla %s evdp %s mag %s \n" % (evlo, evla, evdp, mag)
        s += "wh \n"
    else:
        s += "r %s \n" % filename
        s += "synchronize \n"
        s += "ch o gmt %s %s %s %s %s \n" % (o.year, jday, o.hour, o.minute, o.second)
        s += "ch allt (0 - &1,o&) iztype IO \n"
        s += "ch evlo %s evla %s evdp %s mag %s \n" % (evlo, evla, evdp, mag)
        s += "wh \n"
    s += "q \n"
    p.communicate(s.encode())
    

def rm_response(**kwargs):

##**************************************##

# This function remove instrument response from seismograms directly created
# from SEED or miniSEED, including remean, retrend and taper. All processes are
# executed in sac.

# Created by Meichen Liu on June 10th, 2019 based on sac-manual-v3.6

##**************************************##

##parameters
# filename	the name of file to remove instrument response
# dirname	The directory where SAC files to be processed are saved
# f1-f4		Frequency limits. Lowpass and highpass to surpress low and high
#		frequencies. The four numbers should satisfy f1<f2<f3<f4. f4 
#		should be smaller than Nyquist frequency (1/2 sampling
#		frequency). f3 cannot be too close to f4. The distance between
#		f2 and f3 should be as large as possible. 
# ofile		the appendix of output filename

    import os
    import sys
    import glob
    import subprocess

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    f1 = kwargs.get('f1')
    f2 = kwargs.get('f2')
    f3 = kwargs.get('f3')
    f4 = kwargs.get('f4')
    ofile = kwargs.get('ofile')

    os.chdir('{}'.format(dirname))
    os.putenv("SAC_DISPLAY_COPYRIGHT","O")

    p = subprocess.Popen(['sac'],stdin=subprocess.PIPE)
    s = "wild echo off \n"

    for sacfile in glob.glob('{}'.format(filename)):
        nw,stn,loc,chn,sac,seis_type = sacfile.split('.')
        pz = glob.glob("SAC_PZs_%s_%s_%s_%s_*" % (nw,stn,chn,loc))
        # multi PZ files are not considered
        if len(pz) != 1:
            print("PZ file error for %s" % sacfile)
        else:
            s += "r %s \n" % sacfile
            s += "rmean; rtr; taper \n"
            s += "trans from pol s %s to none freq %s %s %s %s\n" % (pz[0],f1,f2,f3,f4)
            s += "mul 1.0e9 \n"
            s += "w %s.%s.%s.%s.%s.%s \n" % (nw,stn,loc,chn,sac,ofile)

    s += "q \n"
    p.communicate(s.encode())

def rotate(**kwargs):

##**************************************##

# This function rotate NEZ components to RTZ components using sac. To rotate,
# each component of NEZ should exits with same delta. Headers should contain 
# STLA, STLO, EVLA, EVLO, which are necessary for GCP rotation. They should 
# have the same kzdata and kztime. And horizontal components also need to be 
# orthogonal to each other. Rotation would be successful if any of the above is
# not met.

# Created by Meichen Liu on June 10th. 2019 based on sac-manual-v3.6

##**************************************##

##parameters
# dirname	The directory where SAC files are saved.
# filename	the name of SAC file to be rotate
    
    import os
    import sys
    import glob
    import subprocess

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')

    os.putenv("SAC_DISPLAY_COPYRIGHT","O")
    os.chdir('{}'.format(dirname))

    # create sets for each station of NEZ components
    sets = set()
    for sacfile in glob.glob("{}".format(filename)):
        nw, stn, loc, chn, sac, seis_type = sacfile.split('.')
        key = '.'.join([nw, stn, loc, chn[0:2]])
        sets.add(key)

    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off\n"
    for key in sets:
        Z = key + "Z.{}.{}".format(sac,seis_type)
        # if vertical component does not exist, loop to the next station. No
        # rotation would be done on this station.
        if not os.path.exists(Z):
            print("%s: Vertical component missing!" % key)
            continue

        # check if horizontal components exist
        if os.path.exists(key + "E.{}.{}".format(sac,seis_type)) and os.path.exists(key + "N.{}.{}".format(sac,seis_type)):
            E = key + "E.{}.{}".format(sac,seis_type)
            N = key + "N.{}.{}".format(sac,seis_type)
        elif os.path.exists(key + "1.{}.{}".format(sac,seis_type)) and os.path.exists(key + "2.{}.{}".format(sac,seis_type)):
            E = key + "E.{}.{}".format(sac,seis_type)
            N = key + "N.{}.{}".format(sac,seis_type)
        else:
            print("%s: Horizontal components missing!" % key)
            continue

        # check if horizontal components are orthogonal
        Ecmpaz = subprocess.check_output(['saclst','cmpaz','f',E]).decode().split()[1]
        Ncmpaz = subprocess.check_output(['saclst','cmpaz','f',N]).decode().split()[1]
        cmpaz_delta = abs(float(Ecmpaz) - float(Ncmpaz))
        if not (abs(cmpaz_delta-90)<=0.01 or abs(cmpaz_delta-270)<=0.01):
            print("%s: cmpaz1=%s, cmpaz2=%s are not orthogonal!" % (key, Ecmpaz, Ncmpaz))
            continue

        # check B, E, DELTA
        Zb, Ze, Zdelta = subprocess.check_output(['saclst','b','e','delta','f',Z]).decode().split()[1::]
        Eb, Ee, Edelta = subprocess.check_output(['saclst','b','e','delta','f',E]).decode().split()[1::]
        Nb, Ne, Ndelta = subprocess.check_output(['saclst','b','e','delta','f',N]).decode().split()[1::]
        
        if not (float(Zdelta) == float(Edelta) and float(Zdelta) == float(Ndelta)):
            print("%s: delta not equal!" % key)
            continue

        # get the max B and min E to be the data window
        begin = max(float(Zb), float(Eb), float(Nb))
        end = min(float(Ze), float(Ee), float(Ne))

        # output filename
        R, T, Z0 = key + 'R.{}.{}'.format(sac,seis_type), key + 'T.{}.{}'.format(sac,seis_type), key + 'Z.{}.{}'.format(sac,seis_type)

        s += "cut %f %f \n" % (begin, end)
        s += " r %s %s \n" % (E, N)
        s += "rotate to gcp \n"
        s += "w %s %s \n" % (R, T)
        s += "r %s \n" % Z
        s += "w %s \n" % Z0 
    s += "q \n"
    p.communicate(s.encode())

    # delete original files
#    for sacfile in glob.glob("*.BH[NEZ].SAC*"):
#        os.unlink(sacfile)

def myfilter(**kwargs):

##**************************************##

# This function filter seismograms to the desired frequency range using sac.

# Created by Meichen Liu on June 12nd, 2019

##**************************************##

##parameters
# filename	the name of sac file
# dirname	the directory files are saved
# filter_type	including "bp","lp","hp"
# f_low		lower bound frequency
# f_high	higher bound frequency
# ofile		the appendix of output file

    import os
    import subprocess
    import glob
    import sys
    import numpy as np

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    filter_type = kwargs.get('filter_type')
    f_low = kwargs.get('f_low')
    f_high = kwargs.get('f_high')
    ofile = kwargs.get('ofile')

    os.chdir('{}'.format(dirname))
    os.putenv('SAC_DISPLAY_COPYRIGHT','O')

    p = subprocess.Popen(['sac'],stdin=subprocess.PIPE)
    s = "wild echo off\n"

    filelist = glob.glob('{}'.format(filename))
    if len(filelist)>20:
        for i in np.arange(round(len(filelist)/20)):
            s += "r %s\n" % filename
            if filter_type == 'bp':
                s += "bp n 4 p 2 c %f %f \n" % (f_low,f_high)
            elif filter_type == 'lp':
                s += "lp n 4 p 2 c %f \n" % f_low
            elif filter_type == 'hp':
                s += "hp n 4 p 2 c %f \n" % f_high
            s += "w append .%s \n" % ofile
        s += "r %s\n" % filelist[i*20+20::]
        if filter_type == 'bp':
            s += "bp n 4 p 2 c %f %f \n" % (f_low,f_high)
        elif filter_type == 'lp':
            s += "lp n 4 p 2 c %f \n" % f_low
        elif filter_type == 'hp':
            s += "hp n 4 p 2 c %f \n" % f_high
        s += "w append .%s \n" % ofile
    else:
        s += "r %s\n" % filename
        if filter_type == 'bp':
            s += "bp n 4 p 2 c %f %f \n" % (f_low,f_high)
        elif filter_type == 'lp':
            s += "lp n 4 p 2 c %f \n" % f_low
        elif filter_type == 'hp':
            s += "hp n 4 p 2 c %f \n" % f_high
        s += "w append .%s \n" % ofile
      
    s += "q \n"
    p.communicate(s.encode())

def my_resample(**kwargs):

##**************************************##

# This function resample the seismograms to desired frequency.

# Created by Meichen Liu on June 16th, 2019

##**************************************##

##parameters
# filename	the name of sac files
# dirname	the directory seismograms are saved
# delta		resampling rate in seconds
# ofile 	the appendix of output file

    import os
    import glob
    import subprocess
    import obspy

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    delta = kwargs.get('delta')
    freq = 1./float(delta)/2.
    ofile = kwargs.get('ofile')

    os.chdir('{}'.format(dirname))
    os.putenv("SAC_DISPLAY_COPYRIGHT","O")

    subprocess.call(['rm *{}'.format(ofile)],shell=True)
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off\n"

    for sacfile in glob.glob('{}'.format(filename)):
        sac_delta = round(obspy.read('{}'.format(sacfile))[0].stats.sac['delta'],3)
        if sac_delta < float(delta):
            s += "r %s \n" % sacfile
            s += "lp c %f \n" % freq
            s += "interpolate delta %s \n" % delta
            s += "w append .%s \n" % ofile
        else:
            print("%s original sampling rate smaller than %f" % (sacfile,delta))

    s += "q \n"
    p.communicate(s.encode())

def taup_mark(**kwargs):

##**************************************##

# This function save the predicted arrival time of desired phase. 

# Created by Meichen Liu on June 10th, 2019

##**************************************##

##parameters
# filename	the name of sac file
# dirname	the directory seismograms are saved
# phase		name of the phase
# model		the model to calculate travel time.
#		"prem", "iasp91","ak135"
# Tn		the number of Tn to store travel time

    import os
    import glob
    import subprocess
    import numpy as np

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    phase = kwargs.get('phase')
    model = kwargs.get('model')
    Tn = kwargs.get('Tn') 
    
    os.chdir('{}'.format(dirname))

    for sacfile in glob.glob('{}'.format(filename)):
        try:
            subprocess.call(['taup_setsac','-mod',model,'-ph','{}-{}'.format(phase,Tn),'-evdpkm','{}'.format(sacfile)])
        except:
            print('Unable to mark arrival time of {} wave on file {}'.format(phase,sacflie))
            subprocess.call(['rm {}'.format(sacfile)],shell=True)
            continue

def table_mark(**kwargs):

##**************************************##

# This function mark arrival time based on self-made timetable.

##**************************************##

##parameters
# dirname	The directory where sacfiles are saved
# filename	the name of sac files
# ofile		the appendix of output files. False means no appendix
# tablef	the timetable file
# Tn		to save arrival time


    import numpy as np
    import subprocess
    import os
    import glob

    dirname = kwargs.get('dirname')
    filename = kwargs.get('filename')
    ofile = kwargs.get('ofile')
    tablef = kwargs.get('tablef')
    Tn = kwargs.get('Tn')

    os.chdir('{}'.format(dirname))
    os.putenv('SAC_DISPLAY_COPYRIGHT','O')

    subprocess.call(['rm *{}'.format(ofile)],shell=True)
    p = subprocess.Popen(['sac'],stdin=subprocess.PIPE)
    s = "wild echo off\n"

    # read time table: distance(deg) depth(km) phase traveltime ray_param(s/deg) takeoff(deg) incident(deg)
    table = np.genfromtxt(tablef)
    for sacfile in glob.glob('{}'.format(filename)):
        b, e, o, evdp, gcarc = subprocess.check_output(['saclst','b','e','o','evdp','gcarc','f',sacfile]).decode().split()[1::]
        gcarc = round(float(gcarc))
        evdp_lower = np.floor(float(evdp)/10)*10
        evdp_upper = evdp_lower + 10
        gcarc_index = np.where(table[:,0] == gcarc)[0]
        evdp_lower_index = np.where(table[:,1] == evdp_lower)[0]
        lower_index = [x for x in gcarc_index if x in evdp_lower_index]
        tt_lower = table[lower_index,2]
        evdp_upper_index = np.where(table[:,1] == evdp_upper)[0]
        upper_index = [x for x in gcarc_index if x in evdp_upper_index]
        tt_upper = table[upper_index,2]
        if tt_lower != 999999 and tt_upper != 999999:
            tt = (float(evdp) - evdp_lower)/10 * (float(tt_upper) -float(tt_lower)) + float(tt_lower)
            if tt>float(b) and tt<float(e):
                s += "r %s\n" % sacfile
                s += "ch %s %f\n" % (Tn, tt)
                s += "wh\n"
                if ofile != "False":
                    s += "w append %s\n" % ofile

    s += "q \n"
    p.communicate(s.encode())        
    
def seis_type(**kwargs):

##**************************************##

# This function transfer seismograms among displacement, velocity and 
# acceleration types.

# Created by Meichen Liu on Jund 17th, 2019

##**************************************##

##parameters
# filename	the name of sac files
# dirname	the directory where files are saved
# ifile		the type of input file, "dis", "vel" or "acc"
# ofile		the type of output file

    import os
    import subprocess
    import numpy as np
    import glob

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    ifile = kwargs.get('ifile')
    ofile = kwargs.get('ofile')

    os.chdir('{}'.format(dirname))
    os.putenv("SAC_DISPLAY_COPYRIGHT","O")

    p = subprocess.Popen(['sac'],stdin=subprocess.PIPE)
    s = "wild echo off\n"

    filelist = glob.glob('{}'.format(filename))
    if len(filelist)>10:
        for i in np.arange(round(len(filelist)/10)):
            s += "r %s \n" % filelist[i*10:i*10+10]
            if ifile == 'dis':
                if ofile == 'vel':
                    s += "dif\n"
                if ofile == 'acc':
                    s += "dif\n"
                    s += "dif\n"
        
            if ifile == 'vel':
                if ofile == 'dis':
                    s += "dif\n"
                if ofile == 'acc':
                    s += "int\n"
        
            if ifile == 'acc':
                if ofile == 'vel':
                    s += "int\n"
                if ofile == 'dis':
                    s += "int\n"
                    s += "int\n"
            s += "w append .%s\n" % ofile
        s += "r %s \n" % filelist[i*10+10::]
        if ifile == 'dis':
            if ofile == 'vel':
                s += "dif\n"
            if ofile == 'acc':
                s += "dif\n"
                s += "dif\n"
    
        if ifile == 'vel':
            if ofile == 'dis':
                s += "dif\n"
            if ofile == 'acc':
                s += "int\n"
    
        if ifile == 'acc':
            if ofile == 'vel':
                s += "int\n"
            if ofile == 'dis':
                s += "int\n"
                s += "int\n"
        s += "w append .%s\n" % ofile
    else:
        s += "r %s \n" % filename
        if ifile == 'dis':
            if ofile == 'vel':
                s += "dif\n"
            if ofile == 'acc':
                s += "dif\n"
                s += "dif\n"
    
        if ifile == 'vel':
            if ofile == 'dis':
                s += "dif\n"
            if ofile == 'acc':
                s += "int\n"
    
        if ifile == 'acc':
            if ofile == 'vel':
                s += "int\n"
            if ofile == 'dis':
                s += "int\n"
                s += "int\n"
        s += "w append .%s\n" % ofile
        
    s += "q \n"
    p.communicate(s.encode())


def peak_flip(**kwargs):

##**************************************##

# This function find the peak point(maximum or minimum) at the vicinity of
# Tn of a seismogram and save to SAC file header

# Created by Meichen Liu on June 17th, 2019

##**************************************##

##parameters
# filename	the name of sac file
# dirname	the directory seismograms are saved
# T_phase	the name of pick time of phase
# t_before	seconds before Tn of the peak picking window
# t_after	seconds after Tn of the peak picking window
# T_name	the number of T to save peak time
# flip		True or False. If True, then flip the seismogram is peak is 
#		negative
# ofile		the appendix of output file

    import os
    import subprocess
    import numpy as np
    import glob
    import obspy

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    T_phase = kwargs.get('T_phase')
    t_before = kwargs.get('t_before')
    t_after = kwargs.get('t_after')
    T_name = kwargs.get('T_name')
    flip = kwargs.get('flip')
    ofile = kwargs.get('ofile')

    os.chdir('{}'.format(dirname))
    os.putenv("SAC_DISPLAY_COPYRIGHT","O")

    p = subprocess.Popen(['sac'],stdin=subprocess.PIPE)
    s = "wild echo off\n"

    for sacfile in glob.glob('{}'.format(filename)):
        T, delta, b, e = subprocess.check_output(['saclst',T_phase,'delta','b','e','f',sacfile]).decode().split()[1::]
        N_before = int(float(t_before)/float(delta))
        N_after = int(float(t_after)/float(delta))
        N_T = int((float(T)-float(b))/float(delta))

        try:
            tr = obspy.read('{}'.format(sacfile))[0]
            window = tr.data[N_T-N_before:N_T+N_after]
            peak_loc = np.argmax(abs(window)) + N_T - N_before
            peak_t = float(b) + float(peak_loc) * float(delta)
        except:
            print("Failed to find max amp within the window of {}".format(sacfile))
            continue
        
        if flip == 'True' and window[np.argmax(abs(window))]<0:
            s += "r %s \n" % sacfile
            s += "mul -1 \n"
            s += "ch %s %f \n" % (T_name, peak_t)
            s += "wh \n"
            s += "w append %s \n" % ofile
        else:
            s += "r %s \n" % sacfile
            s += "ch %s %f \n" % (T_name, peak_t)
            s += "wh \n"
            s += "w append %s \n" % ofile
#            print("Failed to peak_flip %s" % sacfile)

    s += "q \n"
    p.communicate(s.encode())

        
def quality_control(**kwargs):

##**************************************##

# This function delete low-quality seismograms.

# Created by Meichen Liu on June 11th. 2019

##**************************************##

##parameters
# filename	the name of sac file
# dirname	the directory seismograms are saved
# snr1_t1	Two windows are chosen to calculate the signal-noise ratio, in 
#		which snr1 is before peak S, and snr2 isafter peak S. Here 
#		snr1_t1 is the number of seconds before peak S of the start 
#		time of snr1.
# snr1_t2	number of seconds before peak S of the end time of snr1.
# snr2_t1	number of seconds after peak S of the start time of snr2.
# snr2_t2	number of seconds after peak S of the end time of snr2.
# average_ratio	the critical value that the ratio of each window average to 
#		peak amplitude should be smaller than.
# max_ratio	the critical value that the ratio of each window maximum to the
#		peak amplitude should be smaller than.
# t_before	time in seconds before the peak S to look for the lowest point
# t_after	time in seconds after the peak S to look for the lowest point
# amp_ratio	the critical value that the ratio of minimum to the maximum
#		amplitude should be smaller than.
# ofile		the appendix of output file

    import os
    import glob
    import subprocess
    import sys
    import obspy
    import numpy as np

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    snr1_t1 = kwargs.get('snr1_t1')
    snr1_t2 = kwargs.get('snr1_t2')
    snr2_t1 = kwargs.get('snr2_t1')
    snr2_t2 = kwargs.get('snr2_t2')
    average_ratio = kwargs.get('average_ratio')
    max_ratio = kwargs.get('max_ratio')
    amp_ratio = kwargs.get('amp_ratio')
    t_before = kwargs.get('t_before')
    t_after = kwargs.get('t_after')
    ofile = kwargs.get('ofile')
    

    os.chdir('{}'.format(dirname))
    mk_path('{}/excluded'.format(dirname))
    os.putenv("SAC_DISPLAY_COPYRIGHT","O")

    for sacfile in glob.glob('{}'.format(filename)):
        flag = 0
        peak_t, delta, b = subprocess.check_output(['saclst','t1','delta','b','f',sacfile]).decode().split()[1::]
        peak_loc = int((float(peak_t)-float(b))/float(delta))
        tr = obspy.read('{}'.format(sacfile))[0]
        peak_amp = tr.data[peak_loc]

        # check window snr1
        snr1_N1 = int((float(peak_t)-float(snr1_t1)-float(b))/float(delta))
        snr1_N2 = int((float(peak_t)-float(snr1_t2)-float(b))/float(delta))
        if snr1_N1 <=0 or snr1_N2 >= len(tr.data):
            sys.exit("Windows exceed the range of %s" % sacfile)
        snr1_average = np.sum(abs(tr.data[snr1_N1:snr1_N2]))/(snr1_N2-snr1_N1)
        snr1_max = max(abs(tr.data[snr1_N1:snr1_N2]))
        
        # check window snr2
        snr2_N1 = int((float(peak_t)+float(snr2_t1)-float(b))/float(delta))
        snr2_N2 = int((float(peak_t)+float(snr2_t2)-float(b))/float(delta))
        if snr2_N1 <=0 or snr2_N2 >= len(tr.data):
            sys.exit("Windows exceed the range of %s" % sacfile)
        snr2_average = np.sum(abs(tr.data[snr2_N1:snr2_N2]))/(snr2_N2-snr2_N1)
        snr2_max = max(abs(tr.data[snr2_N1:snr2_N2]))

        # exclude those low-snr seismograms
        if snr1_average > float(average_ratio)*peak_amp or snr2_average > float(average_ratio)*peak_amp:
            flag = 1
        elif snr1_max/peak_amp > float(max_ratio) or snr2_max/peak_amp > float(max_ratio):
            flag = 1

        # exclude those have two wiggles
        N_before = int(float(t_before)/float(delta))
        N_after = int(float(t_after)/float(delta))
        low_window = tr.data[peak_loc-N_before:peak_loc+N_after]
        low_amp = min(low_window)
        if low_amp <0 and (-1)*low_amp/peak_amp >float(amp_ratio):
            flag = 1
#            os.rename(sacfile,'./excluded/{}'.format(sacfile))
        
        # rename those with good quality
        if flag == 0:
            os.rename(sacfile,sacfile+'.{}'.format(ofile))

def cut_seismograms(**kwargs):

##**************************************##

# This function cut seismograms in SAC.

# Created by Meichen Liu on June 17th, 2019

##**************************************##

##parameters
# filename	the name of SAC files
# dirname	the directory seismograms are saved
# tpick		the reference tpick time
# cut_before	seconds cut before tpick
# cut_after	seconds cut after tpick
# ofile		appendix to the output file

    import os
    import glob
    import subprocess
    import numpy as np

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    cut_before = kwargs.get('cut_before')
    cut_after = kwargs.get('cut_after')
    cut_name = kwargs.get('cut_name')
    tpick = kwargs.get('tpick')

    os.chdir('{}'.format(dirname))
    os.putenv("SAC_DISPLAY_COPYRIGHT","O")

    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off\n"
    s += "cut %s -%f %f \n" % (tpick, cut_before, cut_after)

    for sacfile in glob.glob('{}'.format(filename)):
        b, e, t, delta = subprocess.check_output(['saclst','b','e',tpick,'delta','f',sacfile]).decode().split()[1::]
        N_before = int(cut_before/float(delta))
        N_after = int(cut_after/float(delta))

        if float(t)-float(cut_before)<float(b) or float(t)+float(cut_after)>float(e):
            print("Cut window exceeds the range of %s" % sacfile)
        else:
            s += "r %s \n" % sacfile
            s += "w append .%s \n" % ofile

    s += "cut off \n"
    s += "q \n"
    p.communicate(s.encode())

def sac2xh(**kwargs):

##**************************************##

# This function convert ultimate sacfiles to xhfiles and rename them

##**************************************##

##parameters
# filename		the name of sacfiles
# dirname		the path of the directory where sacfiles are 
# seis_type		'dis', 'vel', or 'acc'
# eventid		the id of the event
# chn 			the channel component
# ofile			the name of output xhfile

    import os
    import glob
    import subprocess

    filename = kwargs.get('filename')
    dirname = kwargs.get('dirname')
    seis_type = kwargs.get('seis_type')
    eventid = kwargs.get('eventid')
    chn = kwargs.get('chn')
    ofile = kwargs.get('ofile')

    os.chdir('{}'.format(dirname))
    mk_path('../XH_files')

    f = open('sac2xh.sh','w')
    f.write("#!/bin/bash\n")

    for sacfile in glob.glob('{}'.format(filename)):
        nw, stn, loc, chn = sacfile.split('.')[0:4]
        f.write("sac_2xh {} {}.{}.{}.{}.{}.{}.xh -l 9\n".format(sacfile,eventid,nw,stn,loc,chn,seis_type))
    f.write("cat {}.*.{}.{}.xh > {}.xh \n".format(eventid,chn,seis_type,ofile))

    f.close()
    subprocess.call(['sh sac2xh.sh'],shell=True)
    subprocess.call(['rm sac2xh.sh'],shell=True)

def remove_files(path):
    
    import os
    import subprocess

    os.chdir('{}'.format(path))
    subprocess.call(['rm *.raw *.filtered *.markS *.markP *.markScS *.marksS *.marks410S *.markS410S *.markS660S *.rsm *.*.*.*.dis.xh *.BHR* *.BHZ.* *.BHE* *.BHN* RESP* SAC_PZ* *.SAC*'],shell=True)
    subprocess.call(['rm ../*.mseed'],shell=True)
    subprocess.call(['rm ../*.dataless'],shell=True)
    subprocess.call(['rm ../../stations/*'],shell=True)

def main():
    
    import numpy as np
    import pandas as pd
    import sys

    path = '/home/meichen/Research/SH_TOMO/USArray'
    save_path = '/home/meichen/work1/SH_TOMO/US_TA'

    saveout = sys.stdout
    saveerr = sys.stderr
    f = open('stdout.log','w')
    sys.stderr = f
    sys.stdout = f

    # Read in events info
    events_cat = pd.read_csv('{}/events_cat/events_cat_300_800_2010_2018_Mw6075.txt'.format(path),skipinitialspace=True,header=None,sep=' ')
    events_cat = np.array(events_cat)
    for eventid in events_cat[:,0]:
        folder_path = '{}/event_{}'.format(save_path,eventid)
        index = list(events_cat[:,0]).index(eventid)
        eventtime = events_cat[index,1]
        eventlat = events_cat[index,2]
        eventlon = events_cat[index,3]
        eventdep = events_cat[index,4]
        eventmag = events_cat[index,5]
        print("Start process event_{}".format(eventid))
        # Download event miniseed files
        download_data(event_info=[eventtime,eventlat,eventlon],shape='cir',cenlat=eventlat,cenlon=eventlon,minrad=50,maxrad=120,start=120,end=3600,channel='BH[NEZ]',client='IRIS',save_path=folder_path)
        # Miniseed to SAC
        mseed2sac(save_path=folder_path,jar_path='/home/meichen/bin')
        # Rename SAC files as nw.stn.loc.chn.SAC.raw
        rename(dirname='{}/waveforms/SAC_files'.format(folder_path))
        # Add event info to SAC files
        add_info(filename='*.raw',dirname='{}/waveforms/SAC_files'.format(folder_path),origin=eventtime,evla=eventlat,evlo=eventlon,evdp=eventdep,mag=eventmag)
        # Remove instrument response
        rm_response(filename='*.raw',dirname='{}/waveforms/SAC_files'.format(folder_path),f1=0.004,f2=0.008,f3=1,f4=3,ofile='dis')
        # Rotate NEZ to RTZ
        rotate(filename='*.BH[NEZ].*.dis',dirname='{}/waveforms/SAC_files'.format(folder_path))
#        # Filter to desired frequency band
#        myfilter(filename='*.BHT.SAC.dis',dirname='{}/waveforms/SAC_files'.format(folder_path),filter_type='lp',f_low=0.1,ofile='filtered')
#        myfilter(filename='*.BHR.SAC.dis',dirname='{}/waveforms/SAC_files'.format(folder_path),filter_type='lp',f_low=0.1,ofile='filtered')
#        myfilter(filename='*.BHZ.SAC.dis',dirname='{}/waveforms/SAC_files'.format(folder_path),filter_type='lp',f_low=0.1,ofile='filtered')
        # Use T0 to save arrival time of P and T1 for S wave
        # Also use T9 for s^410S wave
#        # T7 for sS, T8 for ScS waves
#        table_mark(dirname='{}/waveforms/SAC_files'.format(folder_path),filename='*.filtered',ofile=".markS",tablef='/home/meichen/Utils/TauP_calculations/TravelTimeTables/tt.prem.S',Tn='T1')
#        table_mark(dirname='{}/waveforms/SAC_files'.format(folder_path),filename='*.markS',ofile=".markP",tablef='/home/meichen/Utils/TauP_calculations/TravelTimeTables/tt.prem.P',Tn='T0')
#        table_mark(dirname='{}/waveforms/SAC_files'.format(folder_path),filename='*.markP',ofile='.markScS',tablef='/home/meichen/Utils/TauP_calculations/TravelTimeTables/tt.prem.ScS',Tn='T8')
#        table_mark(dirname='{}/waveforms/SAC_files'.format(folder_path),filename='*.markScS',ofile='.marksS',tablef='/home/meichen/Utils/TauP_calculations/TravelTimeTables/tt.prem.sS',Tn='T7')
#        table_mark(dirname='{}/waveforms/SAC_files'.format(folder_path),filename='*.marksS',ofile='.markS410S',tablef='/home/meichen/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.S^410S',Tn='T6')
#        table_mark(dirname='{}/waveforms/SAC_files'.format(folder_path),filename='*.markS410S',ofile='.markS660S',tablef='/home/meichen/Utils/TauP_calculations/TravelTimeTables/tt.prem_ml.S^660S',Tn='T5')
        # Resample seismograms
        my_resample(filename='*.BHT.SAC.dis',dirname='{}/waveforms/SAC_files'.format(folder_path),delta=0.1,ofile='rsm')
        # Create one xh file for one component of one event
#        sac2xh(filename='*.BHZ.*.rsm',dirname='{}/waveforms/SAC_files'.format(folder_path),eventid=eventid,seis_type='dis',chn='BHZ',ofile='{}.BHZ.dis'.format(eventid))
        sac2xh(filename='*.BHT.*.rsm',dirname='{}/waveforms/SAC_files'.format(folder_path),eventid=eventid,seis_type='dis',chn='BHT',ofile='{}.BHT.dis'.format(eventid))
#        sac2xh(filename='*.BHR.*.rsm',dirname='{}/waveforms/SAC_files'.format(folder_path),eventid=eventid,seis_type='dis',chn='BHR',ofile='{}.BHR.dis'.format(eventid))
        # remove unnecessary files
        remove_files('{}/waveforms/SAC_files'.format(folder_path))

    sys.stdout = saveout
    sys.stderr = saveerr
    f.close()

main()
