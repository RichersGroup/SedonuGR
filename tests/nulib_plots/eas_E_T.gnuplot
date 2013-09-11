k_b = 8.6173324e-11  #MeV/K
set grid x
set grid y
set xtics out
set ytics out
filename = "./eas_E_T.dat"
set title "`head -n1 ./eas_E_T.dat`"
set xlabel "Neutrino Energy (MeV)"
set y2label "Colors: log10[Temperature(MeV)]"
set key off

set term unknown
unset log x
unset log y
plot filename using 2:3:(log10($1)) w l palette
set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
set log x
set log y
set ylabel "Integrated Emis (erg/s/cm^3/ster)"
set term pdf
set output "emis_E_T.pdf"
replot
set output

set term unknown
unset log x
unset log y
plot filename using 2:4:(log10($1)) w l palette
set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
set log x
set log y
set ylabel "Absorption Opacity (cm^2/g)"
set term pdf
set output "absopac_E_T.pdf"
replot
set output

set term unknown
unset log x
unset log y
plot filename using 2:5:(log10($1)) w l palette
set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
set log x
set log y
set ylabel "Scattering Opacity (cm^2/g)"
set term pdf
set output "scatopac_E_T.pdf"
replot
set output