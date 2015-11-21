k_b = 1.3806488e-16 # erg/K
c = 2.99e10         #cm/s
h = 6.62606957e-27  #erg/Hz
MeV = 6.774e-6      #erg
sigma = 5.67051e-5  #erg/cm^2/K^4/s
pi = 3.14159265359
visc_specific_heat_rate = 1e22

r = 1e5                            #cm
T = 5e+10                         #K
rho = 1e11                         #g/cm^3
Lvisc = 4./3.*pi*r*r*r*rho * visc_specific_heat_rate #erg/s
print "Lvisc=",Lvisc
Lfermi = 6.*pi*7./8.*sigma*T*T*T*T         #erg/s/cm^2
print "Lfermi=",Lfermi
N=Lvisc/Lfermi*0

set xlabel "Neutrino Frequency (Hz) (2.5e20 Hz/MeV)"
set ylabel "Energy Flux (erg/s/Hz/sr)"
set grid

set xrange [0:5e22]
set term pdf
set output "compare_spectrum_0.pdf"
set title "Electron Neutrinos"
plot N*x*x*x*h/c/c/(exp(h*x/(k_b*T))+1.), './spectrum_species0_00001.dat' using 1:3, '0.dat' u 1:3 w l
set output

set term pdf
set output "compare_spectrum_1.pdf"
set title "Electron Anti-Neutrinos"
plot N*x*x*x*h/c/c/(exp(h*x/(k_b*T))+1.), './spectrum_species1_00001.dat' using 1:3, '1.dat' u 1:3 w l
set output

set term pdf
set output "compare_spectrum_2.pdf"
set title "Mu/Tau Anti/Neutrinos"
plot N*x*x*x*h/c/c/(exp(h*x/(k_b*T))+1.)*4., './spectrum_species2_00001.dat' using 1:3, '2.dat' u 1:3 w l
set output

# planck function has units of erg/s/cm^2/Hz/ster
# 4*pi is the total solid angle neutrinos are collected from
# pi*r*r is the cross-sectional surface area
# 7/8 above in P because Integrate[fermi-dirac function, 0, infinity] = 7/8 * Integrate[planck_function, 0, infinity]
# the 10 in N is because the spectrum was integrated over 10 iterations