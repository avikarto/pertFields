from math import pi, sqrt
import numpy as np  # This is used by other modules importing * from here

# ----------physical constants----------
a0 = 5.2917721092e-11
GEVperHartree = 27.21138505e-9  # GeV per Hartree
eVperHartree = 27.21138505  # eV per Hartree
secPerAU = 2.41888e-17  # seconds per au time
c = 137.035999074  # speed of light
c2, c3 = c**2, c**3
mu0 = 4.0*pi/c2
eps0 = 1.0/(4.0*pi)
eta0 = sqrt(mu0/eps0)  # impedence of free space


# ----------global constant parameters (in a.u. when applicable)----------
mm = 1  # angular L-G index
nn = 1  # radial L-G index.
pertOrder = 1  # = j_max in the BGV sum, not 2j

phi0 = 3.0*pi/2.0  # initial phase
s = 70  # cycles = s-value:   3=64, 5=178, 6=256, 7=349, 8=456, 9=577, 10=712
lambda0 = 800.e-9/a0
w0 = 785.e-9/a0

omega0 = 2.0*pi*c/lambda0  # laser frequency (angular)
k = 2.0*pi/lambda0  # wave number
zr = 0.5*k*w0**2  # Rayleigh range

# polarization choice {1:linear(x), 2:radial}
polar = 2

grid = 200  # number of x,y grid cells for plotting and normalizing.
xmin = -3.0*w0
xmax = -xmin
ymin, ymax = xmin, xmax
dx = (xmax-xmin)/float(grid)
dy = dx

mmReal = float(mm)
nnReal = float(nn)
sReal = float(s)
ii = 1.0j  # imaginary i
epsilonc2 = c/(2.0*zr*omega0)  # \epsilon_c^2
