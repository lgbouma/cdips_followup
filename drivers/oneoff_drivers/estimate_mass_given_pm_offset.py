"""
There are a few cases, like TOI 861.01 and TIC 859480036.01, for which the
proper motion is a little bit off from cluster mean, but far enough that it
seems noticeable & cluster membership is still good.

This suggests motion from the orbital companion, and in turn can be used to
guess the mass.

mu [mas/yr] ~= v_orb [AU/time] / distance [pc]
"""

from astropy import units as u
from math import pi, sqrt
from astropy import constants as const

name = 'TOI861.01'
P = 13.47*u.day
Mtot = 2.66*u.Msun
dist = (1/(2.34e-3))*u.pc
pm_offset = sqrt(2)*u.mas/u.yr

mu_bound = pm_offset
v_orb_bound = mu_bound * dist

K_bound = v_orb_bound.to(u.rad*u.km/u.s)

m2sini_bound = (
    (K_bound.value*u.km/u.s)/(28.4*u.m/u.s) *
    (Mtot.to(u.Msun).value)**(2/3) * (P/(1*u.yr))**(1/3)
)

print('For {}'.format(name))
print('with PM offset of {:.2f}'.format(pm_offset))
print('dist =  {:.2f}'.format(dist))
print('period =  {:.2f}'.format(P))
print('Mtot =  {:.2f}'.format(Mtot))

print('get vorb_bound ~= {:.2f}'.format(v_orb_bound.to(u.rad*u.km/u.s)))
print('get m2sini/mjup ~= {:.2f}'.format(m2sini_bound.cgs))
