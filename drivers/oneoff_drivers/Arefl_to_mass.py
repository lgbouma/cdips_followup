# Eq9 of Shporer 2017, assumes secondary does not give any light

from astropy import units as u, constants as c
import numpy as np

M1 = 8.51*u.Msun
P = 4.6*u.day
R2 = 1.3*u.Rsun
Arefl = 0.015 # i.e. 1.5%
alpha_refl = 1
sini = 1

A_pred = (
    57*alpha_refl*sini
    * (M1/(1*u.Msun))**(-2/3)
    * (P/(1*u.day))**(-4/3)
    * (R2/(1*u.Rjup))**2
)/1e6

print(A_pred.cgs)
