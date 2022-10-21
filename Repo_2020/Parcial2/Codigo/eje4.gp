set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5

set title 'Registro del detector en función del tiempo'
set xlabel 'Tiempo [clicks]'
set ylabel 'Densidad [masa/celda^2]'

set terminal jpeg enhanced
set output 'eje4_detector.jpg'
plot 'eje4_detector.dat' w l ls 1 t 'Detector'

unset title
unset xlabel
unset ylabel
unset grid

set title 'Distribución de densidad en t=5000 clicks'
set xlabel 'ix'
set ylabel 'iy'
set xrange [0:256]

set pm3d map
set size ratio -1
set output 'eje4_5000.jpg'
splot 'eje4_5000.dat'


