set pointsize 0.25
k_b = 1.3806488e-16 # erg/K
c = 2.99e10         #cm/s
h = 6.62606957e-27  #erg/Hz
pi = 3.14159265359
k_MeV = 1.16046e10
MeV = 0.0000016021773

r = 2e5                             #cm
T = 10*k_MeV                        #K
mu = 10*MeV

set xlabel "Neutrino Frequency (Hz) (2.5e20 Hz/MeV)"
set ylabel "Energy Flux (erg/s/Hz/sr)"

set term pdf
set output "compare_spectrum_0.pdf"
set title "Electron Neutrinos"
plot pi*r*r*x*x*x*h/c/c*1/(exp((h*x-mu)/(k_b*T))+1.), './spectrum_species0_00001.dat' using 1:3
set output

set term pdf
set output "compare_spectrum_1.pdf"
set title "Electron Anti-Neutrinos"
plot pi*r*r*x*x*x*h/c/c*1/(exp((h*x+mu)/(k_b*T))+1.), './spectrum_species1_00001.dat' using 1:3
set output

set term pdf
set output "compare_spectrum_2.pdf"
set title "Mu/Tau Anti/Neutrinos"
plot 4.*pi*r*r*x*x*x*h/c/c*1/(exp(h*x/(k_b*T))+1.), './spectrum_species2_00001.dat' using 1:3
set output

# planck function has units of erg/s/cm^2/Hz/ster
# pi*r*r is the cross-sectional surface area
