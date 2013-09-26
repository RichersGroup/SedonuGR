set term pdf
set grid
set xrange [:9.999e5]

set xlabel "Radius (km)"
set ylabel "Y_e"
set output "r_ye.pdf"
plot 'fluid_00000' u 1:7 w l, 'fluid_00001' u 1:7 w l, 'fluid_00002' u 1:7 w l, 'fluid_00003' u 1:7 w l, 'fluid_00004' u 1:7 w l, 'fluid_00005' u 1:7 w l, 'fluid_00006' u 1:7 w l, 'fluid_00007' u 1:7 w l, 'fluid_00008' u 1:7 w l, 'fluid_00009' u 1:7 w l
set output

set xlabel "Radius (km)"
set ylabel "Temperature (K)"
set output "r_T.pdf"
plot 'fluid_00000' u 1:6 w l, 'fluid_00001' u 1:6 w l, 'fluid_00002' u 1:6 w l, 'fluid_00003' u 1:6 w l, 'fluid_00004' u 1:6 w l, 'fluid_00005' u 1:6 w l, 'fluid_00006' u 1:6 w l, 'fluid_00007' u 1:6 w l, 'fluid_00008' u 1:6 w l, 'fluid_00009' u 1:6 w l
set output
