set grid
set style line 1 lc rgb "blue" lt 1 lw 1 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 1 pt 7 ps 0.5
set style line 3 lc rgb "black" lt 1 lw 1 pt 7 ps 0.5

set title 'Trayectoria del electrón en el espacio'
#set tics 0.2
set xlabel 'X [m new]'
set ylabel 'Y [m new]'
set zlabel 'Z '
set view 60,300
set term pdf


set out 'espacio.pdf'
splot 'espacio.dat' w l ls 1 t 'electrón'
