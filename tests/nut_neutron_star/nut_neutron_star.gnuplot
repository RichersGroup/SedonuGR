set term pdf
set grid
set xrange [:9.999e5]

set xlabel "Radius (km)"
set ylabel "Y_e"
set output "r_ye.pdf"
plot 'ray_00000' u 1:4 w l, 'ray_00001' u 1:4 w l, 'ray_00002' u 1:4 w l, 'ray_00003' u 1:4 w l, 'ray_00004' u 1:4 w l, 'ray_00005' u 1:4 w l, 'ray_00006' u 1:4 w l, 'ray_00007' u 1:4 w l, 'ray_00008' u 1:4 w l, 'ray_00009' u 1:4 w l
set output

set xlabel "Radius (km)"
set ylabel "Temperature (K)"
set output "r_T.pdf"
plot 'ray_00000' u 1:3 w l, 'ray_00001' u 1:3 w l, 'ray_00002' u 1:3 w l, 'ray_00003' u 1:3 w l, 'ray_00004' u 1:3 w l, 'ray_00005' u 1:3 w l, 'ray_00006' u 1:3 w l, 'ray_00007' u 1:3 w l, 'ray_00008' u 1:3 w l, 'ray_00009' u 1:3 w l
set output
