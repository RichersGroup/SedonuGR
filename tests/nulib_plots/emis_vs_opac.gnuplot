k_b = 8.6173324e-11  #MeV/K
set grid x
set grid y
set xtics out
set ytics out
filename = "./emis_vs_opac.dat"
set title "`head -n1 ./emis_vs_opac.dat`"
set xlabel "Neutrino Energy (MeV)"
set ylabel "Emis (erg/s/cm^3/ster/Hz)"

set term unknown
plot filename using 1:2 w l title "blackbody*absopacity", filename using 1:3 w l title "emissivity"
set term pdf
set output "emis_vs_abs.pdf"
replot
set output
