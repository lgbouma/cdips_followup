from astropy import units as u, constants as c
import numpy as np

vsini = 167.5*u.km/u.s # SM estimate. pretty wrong!

Prot = 18.6*u.hr
rstar = 0.37*u.Rsun
v = 2*np.pi*rstar / Prot
print( (vsini / v).cgs)
