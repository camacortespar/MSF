set xlabel 'ix'
set ylabel 'iy'
set xrange [0:256]

set title 'Distribución de Densidad Inicial'
set pm3d map
set size ratio -1
set terminal jpeg enhanced
set output 'dens_inicial.jpg'
splot 'p1.dat'