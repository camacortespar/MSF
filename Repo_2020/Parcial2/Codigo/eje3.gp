set xlabel 'ix'
set ylabel 'iy'
set xrange [0:256]

set pm3d map
set size ratio -1
set terminal jpeg enhanced
set output 'eje3_500.jpg'
set title 'Distribución de densidad en t=500 clicks'
splot 'eje3_500.dat'
unset title

set title 'Distribución de densidad en t=1000 clicks'
set output 'eje3_1000.jpg'
splot 'eje3_1000.dat
unset title

set title 'Distribución de densidad en t=2000 clicks'
set output 'eje3_2000.jpg'
splot 'eje3_2000.dat
