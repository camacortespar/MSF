set termoption enhanced
set grid
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 ps 0.5

set title 'Respuesta del Detector'
set xlabel 't (clicks)'
set ylabel 'Densidad'

set terminal jpeg enhanced
set output 'detector.jpg'
plot 'det.dat' w l ls 1 t 'Detector'