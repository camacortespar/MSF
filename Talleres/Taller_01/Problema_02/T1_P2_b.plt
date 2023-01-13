set xlabel "lambda"
set ylabel "R(lambda)"

plot "T1_P2_bv3.dat" using 1:2 title "R vs Lambda"

set term png size 1250, 1000
set output "T1_P2_b.png"
replot
