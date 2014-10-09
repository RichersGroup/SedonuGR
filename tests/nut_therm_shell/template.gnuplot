set log x
unset log y
set title "Solve Ye and T (T=TEMP_HEREMeV Ye=YE_HERE)"
set xlabel "Input Density (g/ccm)"
set ylabel "Output Ye"
set term pdf
set output "both_ye_v_rho.pdf"
plot 'results.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$9:1/0) title "Sedonu", \
     'results.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$3:1/0) w l title "Predicted"

set log x
set log y
set title "Solve Ye and T (T=TEMP_HEREMeV Ye=YE_HERE)"
set xlabel "Input Density (g/ccm)"
set ylabel "Output Temperature (MeV)"
set term pdf
set output "both_T_v_rho.pdf"
plot 'results.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$8:1/0) title "Sedonu", \
     'results.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$2:1/0) w l title "Predicted"

set log x
unset log y
set title "Solve Ye and T (rho=RHO_HEREg/ccm Ye=YE_HERE)"
set xlabel "Input Temperature (MeV)"
set ylabel "Output Ye"
set term pdf
set output "both_ye_v_T.pdf"
plot 'results.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$9:1/0) title "Sedonu", \
     'results.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$3:1/0) w l title "Predicted"

set log x
set log y
set title "Solve Ye and T (rho=RHO_HEREg/ccm Ye=YE_HERE)"
set xlabel "Input Temperature (MeV)"
set ylabel "Output Temperature (MeV)"
set term pdf
set output "both_T_v_T.pdf"
plot 'results.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$8:1/0) title "Sedonu", \
     'results.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$2:1/0) w l title "Predicted"

unset log x
unset log y
set title "Solve Ye and T (rho=RHO_HEREg/ccm T=TEMP_HEREMeV)"
set xlabel "Input Ye"
set ylabel "Output Ye"
set term pdf
set output "both_ye_v_ye.pdf"
plot 'results.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$9:1/0) title "Sedonu", \
     'results.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$3:1/0) w l title "Predicted"

unset log x
set log y
set title "Solve Ye and T (rho=RHO_HEREg/ccm T=TEMP_HEREMeV)"
set xlabel "Input Ye"
set ylabel "Output Temperature (MeV)"
set term pdf
set output "both_T_v_ye.pdf"
plot 'results.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$8:1/0) title "Sedonu", \
     'results.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$2:1/0) w l title "Predicted"


set log x
unset log y
set title "Solve Ye (T=TEMP_HEREMeV Ye=YE_HERE)"
set xlabel "Input Density (g/ccm)"
set ylabel "Output Ye"
set term pdf
set output "ye_v_rho.pdf"
plot 'results_ye.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$9:1/0) title "Sedonu", \
     'results_ye.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$3:1/0) w l title "Predicted"

set log x
set log y
set title "Solve T (T=TEMP_HEREMeV Ye=YE_HERE)"
set xlabel "Input Density (g/ccm)"
set ylabel "Output Temperature (MeV)"
set term pdf
set output "T_v_rho.pdf"
plot 'results_T.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$8:1/0) title "Sedonu", \
     'results_T.dat' u 1:((abs($2-TEMP_HERE)==0)&&(($3-YE_HERE)==0)?$2:1/0) w l title "Predicted"

set log x
unset log y
set title "Solve Ye (rho=RHO_HEREg/ccm Ye=YE_HERE)"
set xlabel "Input Temperature (MeV)"
set ylabel "Output Ye"
set term pdf
set output "ye_v_T.pdf"
plot 'results_ye.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$9:1/0) title "Sedonu", \
     'results_ye.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$3:1/0) w l title "Predicted"

set log x
set log y
set title "Solve T (rho=RHO_HEREg/ccm Ye=YE_HERE)"
set xlabel "Input Temperature (MeV)"
set ylabel "Output Temperature (MeV)"
set term pdf
set output "T_v_T.pdf"
plot 'results_T.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$8:1/0) title "Sedonu", \
     'results_T.dat' u 2:((abs($1-RHO_HERE)==0)&&(($3-YE_HERE)==0)?$2:1/0) w l title "Predicted"

unset log x
unset log y
set title "Solve Ye (rho=RHO_HEREg/ccm T=TEMP_HEREMeV)"
set xlabel "Input Ye"
set ylabel "Output Ye"
set term pdf
set output "ye_v_ye.pdf"
plot 'results_ye.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$9:1/0) title "Sedonu", \
     'results_ye.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$3:1/0) w l title "Predicted"

unset log x
set log y
set title "Solve T (rho=RHO_HEREg/ccm T=TEMP_HEREMeV)"
set xlabel "Input Ye"
set ylabel "Output Temperature (MeV)"
set term pdf
set output "T_v_ye.pdf"
plot 'results_T.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$8:1/0) title "Sedonu", \
     'results_T.dat' u 3:((abs($1-RHO_HERE)==0)&&(($2-TEMP_HERE)==0)?$2:1/0) w l title "Predicted"

