#!/bin/bash

lon=-95
lat=40

gmt begin evtmap eps,png
gmt basemap -Rg -JE$lon/$lat/5i -B
gmt coast -W0.75p -Dc -A10000
gawk '{print $1,$3,$4}' evtinfo_dis | uniq | gawk '{print $3,$2}' | gmt plot -Sa0.3c -W0.01p,black -Ggreen -N
equal_dist -l $lat $lon -g 30 | gmt plot -W0.3p,black,-- -N
equal_dist -l $lat $lon -g 60 | gmt plot -W0.3p,black,-- -N
equal_dist -l $lat $lon -g 90 | gmt plot -W0.3p,black,-- -N
equal_dist -l $lat $lon -g 120 | gmt plot -W0.3p,black,-- -N
equal_dist -l $lat $lon -g 150 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 0 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 30 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 60 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 90 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 120 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 150 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 180 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 210 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 240 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 270 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 300 | gmt plot -W0.3p,black,-- -N
equal_az -l $lat $lon -g 330 | gmt plot -W0.3p,black,-- -N
gmt end
