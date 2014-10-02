set log x
set log y
set title "Varying Density"
set xlabel "Density (g/ccm)"
set ylabel "BB Energy Density (erg/ccm)"
set term pdf
set output "rho.pdf"
plot 'results.dat'   u 2:((abs($3-TEMP_HERE)/$3<1e-4)&&(abs($4-YE_HERE)/$4<1e-4)?$1:1/0) w l title "Sedonu", \
     'predicted.dat' u 2:((abs($3-TEMP_HERE)/$3<1e-4)&&(abs($4-YE_HERE)/$4<1e-4)?$1:1/0) w l title "Predicted"

set log x
set log y
set title "Varying Temperature"
set xlabel "Temperature (MeV)"
set ylabel "BB Energy Density (erg/ccm)"
set term pdf
set output "temp.pdf"
plot 'results.dat'   u 3:((abs($2-RHO_HERE)/$2<1e-4)&&(abs($4-YE_HERE)/$4<1e-4)?$1:1/0) w l title "Sedonu", \
     'predicted.dat' u 3:((abs($2-RHO_HERE)/$2<1e-4)&&(abs($4-YE_HERE)/$4<1e-4)?$1:1/0) w l title "Predicted"

unset log x
set log y
set title "Varying Electron Fraction"
set xlabel "Ye"
set ylabel "BB Energy Density (erg/ccm)"
set term pdf
set output "ye.pdf"
plot 'results.dat'   u 4:((abs($3-TEMP_HERE)/$3<1e-4)&&(abs($2-RHO_HERE)/$2<1e-4)?$1:1/0) w l title "Sedonu", \
     'predicted.dat' u 4:((abs($3-TEMP_HERE)/$3<1e-4)&&(abs($2-RHO_HERE)/$2<1e-4)?$1:1/0) w l title "Predicted"