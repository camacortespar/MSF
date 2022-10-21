set terminal png size 500,500
set size square
set xrange [-1500:1500]
set yrange [-1500:1500]
set grid
set output '../img/SJT_c.png'
set title 'Sol Jupiter Troyano'
plot '../S-J-T_nonperturbated.dat' using 3:4 w l title 'Jupiter', '../S-J-T_nonperturbated.dat' using 5:6 w l title 'Troyano',  '../S-J-T_nonperturbated.dat' using 1:2 w l title 'Sol',
