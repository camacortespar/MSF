set xlabel 'ix'
set ylabel 'iy'
set xrange [0:256]

set pm3d map
set size ratio -1
set terminal jpeg enhanced
set output 'dens_5000.jpg'
set title 'Distribucion de Densidad (t = 5000 clicks)'
splot 'p4_5000.dat'