set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 1 pt 7 ps 0.5

set title 'Varianza contra el tiempo'
set xlabel 'Tiempo [celdas]'
set ylabel 'Varianza'

f(x)=a*x+b
fit f(x) 'eje2.dat' via a,b
set terminal jpeg enhanced
set output 'eje2.jpg'


Fit = sprintf(" {/:Bold Parámetros de regresión} \n Var = a*T+b \n a = %g +/- %g \n b = %g +/- %g", a, a_err, b, b_err)

set obj 2 rect from graph 1, 0 to graph 0.6, 0.25 fc rgb "white"
set lab 2 Fit at graph 0.6, 0.2


plot 'eje2.dat' w l ls 1 t 'Varianza', f(x) ls 2 t'Fit'
