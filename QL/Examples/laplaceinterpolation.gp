# generate data
cd '/home/peter/quantlib/QuantLib/Examples'
# mode 0 = random points, 1 = whole col (center), 2 = whole col (border)
! MODE=2 N=50 DELPERC=0.25 ./laplaceinterpolation >laplaceinterpolation.dat

# 2d heat map
reset
set view map
set pm3d corners2color c1
set cbrange [-1:1]
unset xtics
unset ytics
unset border
unset key
unset colorbox
splot 'laplaceinterpolation.dat' u 1:2:3 w p pt 5 ps 3 palette # original
splot 'laplaceinterpolation.dat' u 1:2:($5==1?-1:$3) w p pt 5 ps 3 palette # destructed
splot 'laplaceinterpolation.dat' u 1:2:4 w p pt 5 ps 3 palette # reconstructued

# 3d plot
reset
set hidden3d
unset key
splot 'laplaceinterpolation.dat' u 1:2:3:3 w p pt 7 ps 0.8 palette # original
splot 'laplaceinterpolation.dat' u 1:2:($3*1/(1.0-$5)):3 w p pt 7 ps 0.8 palette # destructed
splot 'laplaceinterpolation.dat' u 1:2:4:4 w p pt 7 ps 0.8 palette # reconstructued
splot 'laplaceinterpolation.dat' u 1:2:($4-$3) w pm3d # error

