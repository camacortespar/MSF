set xlabel 'ix'
set ylabel 'iy'

set title 'Adveccion - Difusion de Pulso Gaussiano'
set pm3d map
set size ratio 2
set terminal jpeg enhanced size 500, 800
set output 'ad_NOfuentes.jpg'
splot 'ad.dat'