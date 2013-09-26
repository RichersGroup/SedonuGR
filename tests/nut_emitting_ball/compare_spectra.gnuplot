k_b = 1.3806488e-16 # erg/K
c = 2.99e10         #cm/s
h = 6.62606957e-27  #erg/Hz
sigma = 5.67051e-5  #erg/cm^2/K^4/s
pi = 3.14159265359

r = 1e5                            #cm
T = 5.8e10                         #K
rho = 1e10                         #g/cm^3
L = 7./8.*4./3.*pi*r*r*r*rho * 1.0 #erg/s
print L
P = 4*pi*r*r*sigma*T*T*T*T         #erg/s
N = L/P                            # normalize BB spectrum to unit net luminosity, then multiply by expected luminosity

set xlabel "Neutrino Frequency (Hz) (2.5e20 Hz/MeV)"
set ylabel "Energy Flux (erg/s/Hz)"

#set xrange [0:3e15]
set term pdf
set output "compare_spectrum_0.pdf"
set title "Electron Neutrinos"
plot N*4.*pi*pi*r*r*2.*x*x*x*h/c/c*1/(exp(h*x/(k_b*T))+1.)/6., './species0_I1.spec' using 1:2
set output

set term pdf
set output "compare_spectrum_1.pdf"
set title "Electron Anti-Neutrinos"
plot N*4.*pi*pi*r*r*2.*x*x*x*h/c/c*1/(exp(h*x/(k_b*T))+1.)/6., './species1_I1.spec' using 1:2
set output

set term pdf
set output "compare_spectrum_2.pdf"
set title "Mu/Tau Anti/Neutrinos"
plot N*4.*pi*pi*r*r*2.*x*x*x*h/c/c*1/(exp(h*x/(k_b*T))+1.)*2./3., './species2_I1.spec' using 1:2
set output

# planck function has units of erg/s/cm^2/Hz/ster
# 4*pi is the total solid angle neutrinos are collected from
# pi*r*r is the cross-sectional surface area
# 7/8 above in P because Integrate[fermi-dirac function, 0, infinity] = 7/8 * Integrate[planck_function, 0, infinity]