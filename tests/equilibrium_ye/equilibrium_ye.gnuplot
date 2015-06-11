rho_max = 2e15  #g/cm^3
r_max   = 10    #km

set term pdf
set grid
set yrange [0:0.6]

set title "Equilibrium Ye with no neutrino emission (T=5MeV)"
set xlabel "Radius (km)"
set ylabel "Y_e"
set output "r_ye.pdf"
plot 'fluid_00001.dat' using ($1/1e5):5 title "NuLib" #, 'eos.dat' using ($1/1e5):3 w l title "EOS"
set output

set xlabel "Density (g/cm^3)"
set ylabel "Y_e"
set output "rho_ye.pdf"
set xrange [:] reverse
plot "fluid_00001.dat" using 3:5 title "NuLib" #, 'eos.dat' using 2:3 w l title "EOS"
set output
