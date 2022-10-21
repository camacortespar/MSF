r = 6378
R = 6379
set title 'N orbits'
#color definitions
set border lw 1.5
set style line 1 lc rgb '#000000' lt 1 lw 0.5
set style line 2 lc rgb '#c0c0c0' lt 2 lw 0.5
unset key; unset border; set tics scale 0
set format ''
set angles degrees
set xyplane at -1
set view 75,20
set lmargin screen 0
set bmargin screen 0
set rmargin screen 1
set tmargin screen 1
set parametric
set urange[0:360]
set vrange[-90:90]
set isosamples 25
set hidden3d
#since we are using Cartesian coordinates, we don't want this
#set mapping spherical
splot   r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) with lines linestyle 1,	  'world_110m.txt' u (R*cos($1)*cos($2)):(R*sin($1)*cos($2)):(R*sin($2)) w l lw 2 lc rgb 'black', 'OP.dat' u 1:2:3 w l lw 3 lc 5
pause mouse
