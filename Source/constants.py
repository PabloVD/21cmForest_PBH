""""
Physical and cosmological constants
Author: Pablo Villanueva Domingo
Started: October 2018
Last update: March 2021
"""

import numpy as np

#----- Constants (from 21cmFAST) -----#

# Misc
hplanck = 6.62606896e-27    #  erg s
Msun = 1.989e33 # g
G = 6.67259e-8  # cm^3 g^-1 s^-2
m_p = 1.6726231e-24 # g
m_ec2 = 0.5109989e6 # eV
m_pc2 = 0.938272e9 # eV
sigmaT = 6.6524e-25 #cm^2
c = 29979245800.0 #cm/s
kB = 8.6173324e-5 # eV/K    #1.380658e-16 # erg/K
alphafine = 1./(137.035999)
e2 = alphafine*hplanck*c/2./np.pi   # erg * cm
hnu_0 = 13.6    # eV
sigmaH = 7.91e-18   # ion cross section, cm^2
cutoff = 10.    # log cutoff in integrals
num_k = 1000    # num of k values in k integrals, set 1000 for dndlnM plot

# Units
nu_over_ev = 1.60217646e-12/hplanck # =2.417e14, 1 eV = nu_over_ev Hz
year_sec = 3.154e7 # 1yr in seconds
MpcToCm = 3.086e24 # 1 Mpc = MpcToCm cm
eVtoerg = 1.60218e-12 # 1 eV = eVtoerg erg
kelvintoerg = 1.38064878066852e-16 # 1 kelvin = kelvintoerg erg
kelvintoeV = kelvintoerg/eVtoerg   # 1 kelvin = kelvintoeV eV

# Cosmo
G_N = 6.67259e-8 # cm^3 g^-1 s^-2
hlittle = 0.678 #0.678
H0 = hlittle*3.2407e-18 # s^-1
Omega_b = 0.02226/hlittle**2.
Omega_m = 0.308 #0.3158
Omega_dm = Omega_m - Omega_b
Omega_l = 0.68 #1.-Omega_m
Y_He = 0.245    # mass fraction
mu = 1./(1.-3.*Y_He/4.) # mean molecular weight, mu=1.22 for neutral H + He
rho_c = 3.*H0**2./(8.*np.pi*G) # g cm^-3
n0_h = rho_c*Omega_b*(1.-Y_He)/m_p   # #/cm^3
n0_he = rho_c*Omega_b*Y_He/(4.*m_p)  # #/cm^3
n0_b = n0_h + n0_he     # n0_b = rho_c/(mu*m_p)*Omega_b
Tcmb0 = 2.7255 # K
rho_cmb0 = 8.*np.pi**5.*(kB*Tcmb0)**4./15./(hplanck*c/eVtoerg)**3. # eV/cm^3
LEdd = 4.*np.pi*G_N*Msun*c*m_p/sigmaT   # LEdd/(M/Msun) in erg/s (7.84635981059e+46 keV/s)
zdec = 150.*(Omega_b*hlittle**2./0.022)**(2./5.) - 1.   # z decoupling of Tk


# 21 cm and Ly alpha
Tstar = 0.0682   # h\nu_21 [K]
A10 = 2.85e-15  # hyperfine Einstein coefficient s^-1
nu21 = 1420.4e6    # hyperfine frecuency Hz
f12 = 0.416     # f12 for Lymanalpha
nualpha = 2.47e15   # Hz
Aalpha = 6.25e8 # s^-1
tauGP = 3.*c**3.*Aalpha*n0_h/(8.*np.pi*nualpha**3.*H0*np.sqrt(Omega_m))*10.**(3./2.)    # tau_GP at z=10 and x_HI=1, da mal!
#print "{:.1e}".format(tauGP), "{:.1e}".format(2.7e4*10.**(3./2.))


#----- Precision Constants -----#

nummass = 50
numtau = 50
numalf = 300

"""
nummass = 20
numtau = 25
numalf = 100
"""
