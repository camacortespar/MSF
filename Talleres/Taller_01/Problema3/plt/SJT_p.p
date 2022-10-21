set terminal png size 500,500
set grid
set output '../img/SJT_d.png'
set title 'Sol Jupiter Troyano-perturbado'
plot "../S-J-T_perturbated.dat" using 3:4 w l title 'Jupiter', "../S-J-T_perturbated.dat" using 5:6 w l title 'Troyano', "../S-J-T_perturbated.dat" using 1:2 w l title 'Sol'
