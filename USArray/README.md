This project use SH reverberations to image mantle layering beneath USArray.

data_download.py:
    Utilize obspy.mass_downloader to download a large amount of data.
    Downloaded data is in miniseed files without event info.
    This script include converting miniseed to SAC, adding event info, remove instrument response, rotating.

xh_process.py:
    Write processes to a bash file. A script to illustrate all processes and parameters.

xh_readfile_distance_recordsection.py:
    Read in a xh file and plot the recordsection of distance VS. time. Seismograms are averaged in a 0.1 degree distance bin. Red and blue represent positive and negative polarities, respectively. Setting the upper bound and lower bound is enabled.
