#!/bin/bash

gmt begin model eps,png
gmt basemap -R-133/-63/22/53 -Jx0.2d -Bxa10 -Bya5
gmt plot -Sc0.05c,black -Gblack grid.txt
gmt plot -Sc0.8c -Wred <<EOF
-110 40
EOF
gmt coast -Wthinnest -Df -A10000 -N1
gmt end

