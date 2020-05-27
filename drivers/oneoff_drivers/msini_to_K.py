"""
Eq 14 of Lovis & Fischer 2010, assuming Mstar >> Mplanet
"""
import numpy as np
from astropy import units as u, constants as const

msini = 1*u.Mjup # 49 * u.Mearth
Mstar = 1.1 * u.Msun
P = 8.3 * u.day
e = 0

K = (
    msini.to(u.Mjup).value * (28.4329*(u.m/u.s)).value
    * (1-e**2)**(-1/2) *
    (Mstar.to(u.Msun).value)**(-2/3) *
    (P.to(u.yr).value)**(-1/3)
)

print('K = {:.4f} m/s'.format(K))
