set xlabel 'ix'
set ylabel 'iy'
set xrange [0:256]

set pm3d map
set size ratio -1
set terminal jpeg enhanced
set output 'dens_500.jpg'
set title 'Distribucion de Densidad (t = 500 clicks)'
splot 'p3_500.dat'
unset title

set title 'Distribucion de Densidad (t = 1000 clicks)'
set output 'dens_1000.jpg'
splot 'p3_1000.dat'
unset title

set title 'Distribucion de Densidad (t = 2000 clicks)'
set output 'dens_2000.jpg'
splot 'p3_2000.dat'