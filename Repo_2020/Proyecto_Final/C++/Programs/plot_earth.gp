set terminal pngcairo
set output 'fig.png'

set xr [-2:2]
set yr [-2:2]
set zr [-2:2]

#color definitions
set border lw 1.5
set style line 1 lc rgb '#000000' lt 1 lw 2
set style line 2 lc rgb '#c0c0c0' lt 2 lw 1

unset key; unset border; set tics scale 0

set format ''
set angles degrees

set xyplane at -1
set view 56,81

set lmargin screen 0
set bmargin screen 0 
set rmargin screen 1
set tmargin screen 1 

set parametric
set isosamples 25
set urange[0:360]
set vrange[-90:90]
r = 0.99
R = 1.00

set hidden3d

#since we are using Cartesian coordinates, we don't want this
#set mapping spherical

splot \
  r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) with lines linestyle 2, \
  'world_110m.txt' u (R*cos($1)*cos($2)):(R*sin($1)*cos($2)):(R*sin($2)) w l lw 2 lc rgb 'black', \
  'OP.dat' u 1:2:3 w l lw 3 lc rgb 'red'