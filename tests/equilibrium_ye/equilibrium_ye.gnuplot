rho_max = 2e15  #g/cm^3
r_max   = 10    #km

set term pdf
set grid

set title "Equilibrium Ye with no neutrino emission (T=5MeV)"
set xlabel "Radius (km)"
set ylabel "Y_e"
set output "r_ye.pdf"
plot 'fluid_00001' using ($1/1e5):7
set output

set xlabel "Density (g/cm^3)"
set ylabel "Y_e"
set output "rho_ye.pdf"
set xrange [:] reverse
plot "fluid_00001" using 4:7
set output
