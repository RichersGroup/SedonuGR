rho_max = 2e15  #g/cm^3
r_max   = 10    #km

set term pdf
set grid

set title "Equilibrium Ye with no neutrino emission"
set xlabel "Radius (km)"
set ylabel "Y_e"
set output "r_ye.pdf"
plot 'profile.dat' using ($1/1e5):4
set output

unset title
set xlabel "Density (g/cm^3)"
set ylabel "Y_e"
set output "rho_ye.pdf"
set xrange [:] reverse
plot "profile.dat" using 6:4
set output