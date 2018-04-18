set pointsize 0.25

set log x
unset log y
set title "(rho=RHO_HEREg/ccm Ye=YE_HERE)"
set xlabel "Input Temperature (MeV)"
set ylabel "Output Ye"
set yrange [0:0.55]
set term pdf
set output "ye_v_T.pdf"
plot 'results.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$3:1/0) w l title "Predicted", \
     'results.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$6:1/0) title "Sedonu"

set log x
set log y
set title "(rho=RHO_HEREg/ccm Ye=YE_HERE)"
set xlabel "Input Temperature (MeV)"
set ylabel "Output Temperature (MeV)"
set yrange [0.1:100]
set term pdf
set output "T_v_T.pdf"
plot 'results.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$2:1/0) w l title "Predicted", \
     'results.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$5:1/0) title "Sedonu"

unset log x
unset log y
set title "(rho=RHO_HEREg/ccm T=TEMP_HEREMeV)"
set xlabel "Input Ye"
set ylabel "Output Ye"
set yrange [0:0.55]
set term pdf
set output "ye_v_ye.pdf"
plot 'results.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$3:1/0) w l title "Predicted", \
     'results.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$6:1/0) title "Sedonu"

unset log x
set log y
set title "(rho=RHO_HEREg/ccm T=TEMP_HEREMeV)"
set xlabel "Input Ye"
set ylabel "Output Temperature (MeV)"
set yrange [0.1:100]
set term pdf
set output "T_v_ye.pdf"
plot 'results.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$2:1/0) w l title "Predicted", \
     'results.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$5:1/0) title "Sedonu"

set log x
unset log y
set title "(T=TEMP_HEREMeV Ye=YE_HERE)"
set xlabel "Input Density (g/ccm)"
set ylabel "Output Ye"
set yrange [0:0.55]
set term pdf
set output "ye_v_rho.pdf"
plot 'results.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$3:1/0) w l title "Predicted", \
     'results.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$6:1/0) title "Sedonu"

set log x
set log y
set title "(T=TEMP_HEREMeV Ye=YE_HERE)"
set xlabel "Input Density (g/ccm)"
set ylabel "Output Temperature (MeV)"
set yrange [0.1:100]
set term pdf
set output "T_v_rho.pdf"
plot 'results.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$2:1/0) w l title "Predicted", \
     'results.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$5:1/0) title "Sedonu"
