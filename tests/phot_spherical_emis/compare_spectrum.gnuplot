k_b = 1.3806488e-16 # erg/K
c = 2.99e10         #cm/s
h = 6.62606957e-27  #erg/Hz
sigma = 5.67051e-5  #erg/cm^2/K^4/s
pi = 3.14159265359

r = 2e15                   #cm
T = 1e4                    #K
L = 1e43                   #erg/s
P = 4*pi*r*r*sigma*T*T*T*T #erg/s
N = L/P

set xlabel "Frequency (Hz)"
set ylabel "Energy Flux (erg/s/Hz)"

set xrange [0:2e15]
set yrange [0:2e28]
set term pdf
set output "compare_spectrum.pdf"
plot N*4*pi*pi*r*r*2*x*x*x*h/c/c*1/(exp(h*x/(k_b*T))-1), './species0_I2.spec' using 1:2
set output

# planck function has units of erg/s/cm^2/Hz/ster
# 4*pi is the total solid angle through which photons are collected
# pi*r*r is the cross-sectional surface area