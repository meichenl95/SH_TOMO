#!/bin/bash

lon=-95
lat=40

gmt begin equal_az pdf,png
gmt basemap -Rg -JE$lon/$lat/5i -B
gmt coast -Wthinnest -Dc -A10000
gawk 'NR<2926{print $4,$3}' events6070.txt | gmt plot -Sa0.2c -W0.01p,black -Gred -N
equal_dist -l $lat $lon -g 30 | gmt plot -W0.01p,brown -N
equal_dist -l $lat $lon -g 60 | gmt plot -W0.01p,brown -N
equal_dist -l $lat $lon -g 90 | gmt plot -W0.01p,brown -N
equal_dist -l $lat $lon -g 120 | gmt plot -W0.01p,brown -N
equal_dist -l $lat $lon -g 150 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 0 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 30 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 60 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 90 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 120 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 150 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 180 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 210 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 240 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 270 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 300 | gmt plot -W0.01p,brown -N
equal_az -l $lat $lon -g 330 | gmt plot -W0.01p,brown -N
gmt end
