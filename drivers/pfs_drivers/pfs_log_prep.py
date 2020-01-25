import numpy as np, pandas as pd
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.utils import (
    get_cdips_candidates,
    given_sourceid_get_gaiarow
)
from numpy import array as nparr
from astropy.coordinates import ICRS
from astropy import units as u

from astroquery.gaia import Gaia

source_id = '5510676828723793920'

gaia_r = given_sourceid_get_gaiarow(source_id)
ra, dec = float(gaia_r['ra']), float(gaia_r['dec'])
pmra, pmdec = float(gaia_r['pmra']), float(gaia_r['pmdec'])

c = ICRS(nparr(ra)*u.deg, nparr(dec)*u.deg)

rahmsstr = c.ra.to_string(u.hour, sep=' ', pad=True)
decdmsstr = c.dec.to_string(u.degree, sep=' ', pad=True)

#Please note that proper motion of RA is in seconds of time per year, as required by Magellan telescope control software, which means: divide the proper motion in arcsec by 15, then divide by cosine of Declination.
# Columns are pipe-separated: Name | RA | Dec | Equinox | pmRA (sec/year) | pmDC (arcsec/year) | Epoch | Binning (1x2 or 3x3) | Speed (Normally “Normal”) | Comment including at least magnitude
pmDC = pmdec / 1e3
pmRA = (pmra/(1e3*15))/np.cos(np.deg2rad(dec))

print('pmRA (sec/year) | pmDC (arcsec/year) | Epoch')
print('{:.5f} | {:.4f} | {:.1f} '.format(
    pmRA, pmDC, 2015.5
))
