set xlabel 'ix'
set ylabel 'iy'
set xrange [0:256]

set title 'Distribución de densidad en t=1 clicks'
set pm3d map
set size ratio -1
set terminal jpeg enhanced
set output 'eje1.jpg'
splot 'eje1.dat'
