set xlabel "r"
set ylabel "R(r)"

plot "T1_P2_d_zero1.dat" using 1:2 title "Lambda = 2.40483", \
     "T1_P2_d_zero2.dat" using 1:2 title "Lambda = 5.52008", \
     "T1_P2_d_zero3.dat" using 1:2 title "Lambda = 8.65373", \
     "T1_P2_d_zero4.dat" using 1:2 title "Lambda = 11.7915", \
     "T1_P2_d_zero5.dat" using 1:2 title "Lambda = 14.9309", \

set term png size 1250, 1000
set output "T1_P2_d.png"
replot

