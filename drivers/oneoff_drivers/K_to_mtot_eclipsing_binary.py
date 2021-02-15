"""
https://en.wikipedia.org/wiki/Binary_mass_function

f = M2^3 sin^3 i / (M1 + M2)^2 = Porb K^3 / (2pi * G)

so, you only get M2 or M1 in limiting extreme cases.

i=90 helps.
"""

import numpy as np
from astropy import units as u, constants as const

K =  87942.55 * u.m/u.s
# M1 = 0.45 * u.Msun #Bp-Rp ~= 2.1
P = 1.52 * u.day
e = 0
sini = 1 # assume

f = P * K**3 / (2*np.pi * const.G)

# f = M2^3 / (M1+M2)^2
#
# assume M1=M2, then
#
# f = M2^3 / (2 M2)^2 = 0.25 * M2

print(f'f = M2^3 / (M1+M2)^2 = {f.to(u.Msun):.4f}')
print(f'If M1 = M2... then M2 = {4*f.to(u.Msun):.2f}')

