set xlabel "r"
set ylabel "R(r)"

plot "T1_P2_c_zero1.dat" using 1:2 title "Lambda = 2.4", \
     "T1_P2_c_zero2.dat" using 1:2 title "Lambda = 5.52", \
     "T1_P2_c_zero3.dat" using 1:2 title "Lambda = 8.65", \
     "T1_P2_c_zero4.dat" using 1:2 title "Lambda = 11.8", \
     "T1_P2_c_zero5.dat" using 1:2 title "Lambda = 14.94", \

set term png size 1250, 1000
set output "T1_P2_c.png"
replot

