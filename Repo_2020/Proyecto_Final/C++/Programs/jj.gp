opaque


set ylabel 'degres'
unset xlabel
plot 'coesOP.dat' u 7:4 w l t'ta'
unset ylabel
set ylabel 'km'
set xlabel 'time(s)'
plot 'coesOP.dat' u 7:1 w l t'a'

unset ylabel
unset xlabel
set format y "10^{%L}"	
set ytics '0.0001'
plot 'coesOP.dat' u 7:2 w l t'e'
unset ytics
set ylabel 'degres'
set xlabel 'time(s)'
plot 'coesOP.dat' u 7:3 w l t'i'
unset ylabel

unset xlabel
plot 'coesOP.dat' u 7:5 w l t'aop'
set xlabel 'time(s)'
plot 'coesOP.dat' u 7:6 w l t'raan'

unset multiplot
pause mouse