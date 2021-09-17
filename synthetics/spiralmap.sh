#!/bin/bash

lat=40
lon=-95

gmt begin spiral png,eps
gmt basemap -Rg -JE$lon/$lat/120/5i -B30
gmt coast -Wthinnest -Dc -A10000
gawk 'NR>1&&NR<14{print $3,$2}' eqs_list.txt | gmt plot -Sa0.3c -W0.05,black -Ggreen -N
gawk 'NR>13{print $3,$2}' eqs_list.txt | gmt plot -Sa0.3c -Gblack -N
#gawk 'NR>1{print $2,$1}' stn_list_1x1.txt | gmt plot -St0.05c -Ggreen -N
gawk 'NR>1{print $2,$1}' stn_list.txt | gmt plot -Sc0.05c -W0.01,black -Gblack -N
#echo $lon $lat | gmt plot -Ss0.15c -Ggreen
equal_dist -l $lat $lon -g 30 | gmt plot -W0.01p,grey -N
equal_dist -l $lat $lon -g 60 | gmt plot -W0.01p,grey -N
equal_dist -l $lat $lon -g 90 | gmt plot -W0.01p,grey -N
equal_az -l $lat $lon -g 0 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 30 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 60 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 90 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 120 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 150 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 180 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 210 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 240 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 270 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 300 | gmt plot -W0.01p,gray -N
equal_az -l $lat $lon -g 330 | gmt plot -W0.01p,gray -N
gmt end
