set terminal png size 500,500
set grid
# Punto a
set output '../img/SJ_a.png'
set multiplot
set title 'Sol-Jupiter'
plot "../S-J_original_axis.dat" using 3:4 title "Jupiter" w l, "../S-J_original_axis.dat" using 1:2 title "Sol" w l
set origin .25,.25
set size square .5,.5
clear
unset key
unset object
set xtics -1,.5,1
set ytics -1,.5,1
plot "../S-J_original_axis.dat" using 1:2 title "Sol" w l linecolor "dark-green"
unset multiplot



