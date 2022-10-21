set terminal png size 1000,500
set grid
# Punto b
set output "../img/SJ_b.png"
set multiplot
set title 'Sol-Jupiter rotado'
plot "../S-J_rotated_axis.dat" using 3:4 w l title "Jupiter"
#, "../S-J_rotated_axis.dat" using 1:2 title "Sol" w l
set origin .15,.5
set size .35,.35
clear
unset key
unset object
unset xtics
unset title
set label "Sol"
plot "../S-J_rotated_axis.dat" using 1:2 w l title "Sol" linecolor "green"
unset multiplot
unset output
unset object
unset origin
unset size
unset title
unset xtics

