
set termoption enhanced
set grid
set size square
set style line 1 lc rgb "blue" lt 1 lw 1 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 1 pt 7 ps 0.5
set style line 3 lc rgb "black" lt 1 lw 1 pt 7 ps 0.5

set title 'Trayectoria del electrón en el plano de fase'
set xlabel 'x [m new]'
set ylabel 'Px [m new]'
set term pdf

set out 'planox.pdf'
plot 'plano.dat' u 1:4 w l ls 1 t 'electrón'

set title 'Trayectoria del electrón en el plano de fase'
set xlabel 'y [m new]'
set ylabel 'Py [m new]'
set term pdf

set out 'planoy.pdf'
plot 'plano.dat' u 2:5 w l ls 1 t 'electrón'

set title 'Trayectoria del electrón en el plano de fase'
set xlabel 'x [m new]'
set ylabel 'Px [m new]'
set term pdf

set out 'planoz.pdf'
plot 'plano.dat' u 3:6 w l ls 1 t 'electrón'