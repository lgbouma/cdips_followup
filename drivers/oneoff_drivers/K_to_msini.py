"""
Eq 14 of Lovis & Fischer 2010, assuming Mstar >> Mplanet
"""
import numpy as np
from astropy import units as u, constants as const

K = 4.4 * u.km/u.s
Mstar = 1.22 * u.Msun
P = 3.4 * u.day
e = 0

msini_mjup = (
    (K / (28.4329*(u.m/u.s))).cgs.value * (1-e**2)**(1/2) *
    (Mstar.to(u.Msun).value)**(2/3) *
    (P.to(u.yr).value)**(1/3)
)

print('msini = {:.2f} mjup'.format(msini_mjup))
