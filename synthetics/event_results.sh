#!/bin/bash

gmt begin event_results pdf,png
gmt set COLOR_FOREGROUND RED4
gmt set COLOR_BACKGROUND BLUE4
gmt basemap -Rg -JH8i -Bxa60 -Bya30
gmt coast -Wthinnest -Df -A10000 -Ggrey -Sgrey
gmt makecpt -Cpolar -T-0.05/0.05 -A50 -Z -H > icpt.cpt
#gmt makecpt -Chaxby -T-0.05/0.05 -A50 -Z -H > icpt.cpt
gawk '{print $3,$2,$4}' evtout.txt | gmt plot -Sc0.2c -Cicpt.cpt
gmt colorbar -Cicpt.cpt -Ba0.05 -Dn0.5/-0.1+jCM+w15c/0.3c+e+h
#xh_stationinfo ~/work1/SH_TOMO/Combined_xh/finalxhfiles/all_all.BHT.vel.snrS.xh | gawk '{print $8,$7}' | gmt plot -St0.05c -Gblack@90
gmt end

gmt begin gridpoint_results pdf,png
gmt set COLOR_FOREGROUND RED4
gmt set COLOR_BACKGROUND BLUE4
gmt basemap -R-130/-65/25/50 -JM5i -Bxa10 -Bya5
gmt xyz2grd gridpointout_0400.txt -Ggridpointout.nc -I1 -R-130/-65/25/50
gmt makecpt -Cpolar -T-0.1/0.1 -A50 -Z -H > icpt.cpt
gmt grdimage gridpointout.nc -Cicpt.cpt
gmt coast -Wthinnest -Df -A10000 -N1 -N2
gmt colorbar -Cicpt.cpt -Ba0.05 -Dn0.5/-0.15+jCM+w12c/0.3c+e+h
gmt end

rm icpt.cpt
rm gridpointout.nc
