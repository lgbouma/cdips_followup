"""
Following Wong 2020, Section5.3, Eqs 14 through 16
"""

from astropy import units as u
import numpy as np

###############################################
# Define quantities used in Eqs 14 through 16 #
###############################################

Ra = 4.92*u.Rsun
Rb = 1.29*u.Rsun
a = 24.42*u.Rsun

Teff_a = (10**4.33)*u.K
Teff_b = (10**3.84)*u.K

beta_a = 10**(
    17 - 4*np.log10(Teff_a.cgs.value) - ((11600*u.K)/Teff_a).cgs.value
)
beta_b = 10**(
    17 - 4*np.log10(Teff_b.cgs.value) - ((11600*u.K)/Teff_b).cgs.value
)

###############################################
# Define quantities used in Eqs 14 through 16 #
###############################################

B_1_ill = (
    17/16 * (Ra/a)**2 * (0.25*Ra/a + 1/3)
    * (Teff_b/Teff_a)**4 * (Rb/Ra)**2
    -
    17/16 * (Rb/a)**2 * (0.25*Rb/a + 1/3)*beta_b/beta_a
)

B_2_ill = (
    17/16 * (Ra/a)**2 * (3/16*Ra/a + 16/(27*np.pi**2))
    * (Teff_b/Teff_a)**4 * (Rb/Ra)**2
    -
    17/16 * (Rb/a)**2 * (3/16*Rb/a + 16/(27*np.pi**2))*beta_b/beta_a

)

print(f'B_1_ill = {B_1_ill.cgs:.2e}')
print(f'B_2_ill = {B_2_ill.cgs:.2e}')
print(f'Bolometric correction term a: {beta_a:.3e}')
print(f'Bolometric correction term b: {beta_b:.3e}')
