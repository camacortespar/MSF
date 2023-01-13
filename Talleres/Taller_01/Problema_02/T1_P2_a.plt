set xlabel "r"
set ylabel "R(r)"

plot "T1_P2_av2.dat" using 1:2 title "R vs r" 
set term png size 1250, 1000
set output 'T1_P2_a.png'
replot
