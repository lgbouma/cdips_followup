import numpy as np
from astropy import units as u

# name = 'TIC29786532.01'
# rp_rs = 0.4704  # measured Rp/Rstar
# Tdur = 2.65*u.hr
# P = 1.588*u.day
# rstar = 0.85*u.Rsun # TIC8
# mstar = 0.53*u.Msun # TIC8

name = 'PTFO8*'
rp_rs = np.sqrt(1.5e-2)  # measured Rp/Rstar
P = 0.4485*u.day
Tdur = (9e-2 * P)
rstar = 1.39*u.Rsun # TIC8
mstar = 0.39*u.Msun # TIC8

name = '390252502'
rp_rs = 0.114  # measured Rp/Rstar
P =  5.437*u.day
Tdur = (0.05 * P)
# Teff 6500K... F5V
rstar = 1.46*u.Rsun # mamajek + Teff + HRD
mstar = 1.33*u.Msun # mamajek + Teff + HRD


rp = (rp_rs * rstar).to(u.Rjup)
rho_rhosun = (Tdur/(13*u.hr))**(-3)* (P/(1*u.yr)) # b = 0


print('From transit lightcurve for {}'.format(name))
print('Tdur = {:.2f}'.format(Tdur.to(u.hr)))
print('Rp/R* = {:.4f}'.format(rp_rs))
print('R* = {} (TIC8)'.format(rstar))
print('Rp = {:.2f} (TIC8 Rstar +transit lightcurve depth)'.format(rp.to(u.Rjup)))
print('(density/solar density) = {:.3f} (Tdur+P from LC, assume b=0)'.
      format(rho_rhosun.to(u.d/u.d)))

print('M* = {} (TIC8)'.format(mstar))
print('(density/solar density) = {:.3f} (TIC8)'.format(
    (mstar/(1*u.Msun) * (rstar/(1*u.Rsun))**(-3))
))
