
reset

# png
set terminal pngcairo size 350,262 enhanced font 'Verdana,10'

# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
system('rm -r png')
unset key
set border 0
set view 342,0
set xrange [-1.5:1.5]
set yrange [-1.5:1.5]
set zrange [-0.5:0.5]
system('mkdir -p png')
# spiral upwards
n=0
do for [ii=1:299] {
    n=n+1
    set output sprintf('png/atomo%03.0f.png',n)
    splot 'espacio.dat' every ::1::ii w l ls 1, \
          'espacio.dat' every ::ii::ii w p ls 1
}
