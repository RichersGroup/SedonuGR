set term pdf
set grid
set xrange [:9.999]

set xlabel "Radius (km)"
set ylabel "Theta"
set zlabel "Y_e"
set output "r_ye.pdf"
splot 'fluid_00000.dat' u ($1/100000):2:6 w l, 'fluid_00001.dat' u ($1/100000):2:6 w l, 'fluid_00002.dat' u ($1/100000):2:6 w l, 'fluid_00003.dat' u ($1/100000):2:6 w l, 'fluid_00004.dat' u ($1/100000):2:6 w l, 'fluid_00005.dat' u ($1/100000):2:6 w l, 'fluid_00010.dat' u ($1/100000):2:6 w l, 'fluid_00015.dat' u ($1/100000):2:6 w l, 'fluid_00020.dat' u ($1/100000):2:6 w l
set output

set xlabel "Radius (km)"
set ylabel "Theta"
set zlabel "Temperature (MeV)"
set output "r_T.pdf"
splot 'fluid_00000.dat' u ($1/100000):2:5 w l, 'fluid_00001.dat' u ($1/100000):2:5 w l, 'fluid_00002.dat' u ($1/100000):2:5 w l, 'fluid_00003.dat' u ($1/100000):2:5 w l, 'fluid_00004.dat' u ($1/100000):2:5 w l, 'fluid_00005.dat' u ($1/100000):2:5 w l, 'fluid_00010.dat' u ($1/100000):2:5 w l, 'fluid_00015.dat' u ($1/100000):2:5 w l, 'fluid_00020.dat' u ($1/100000):2:5 w l
set output
