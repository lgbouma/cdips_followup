import numpy as np, pandas as pd
from cdips_followup.utils import (
    given_sourceid_get_gaiarow
)
from numpy import array as nparr
from astropy.coordinates import ICRS
from astropy import units as u

from astroquery.gaia import Gaia

source_id = '2103737241426734336'

def get_log_row(source_id, whichgaia='gaiaedr3'):

    gaia_r = given_sourceid_get_gaiarow(source_id, whichgaia=whichgaia)
    ra, dec = float(gaia_r['ra']), float(gaia_r['dec'])
    pmra, pmdec = float(gaia_r['pmra']), float(gaia_r['pmdec'])

    c = ICRS(nparr(ra)*u.deg, nparr(dec)*u.deg)

    rahmsstr = c.ra.to_string(u.hour, sep=' ', pad=True)
    decdmsstr = c.dec.to_string(u.degree, sep=' ', pad=True)
    signstr = '+' if c.dec > 0 else '-'
    decdmsstr = signstr+decdmsstr

    G = float(gaia_r['phot_g_mean_mag'])
    Bp = float(gaia_r['phot_bp_mean_mag'])
    Rp = float(gaia_r['phot_rp_mean_mag'])

    # Eq 3 Stassun+2019, which is from Evans et al., 2018. DR2.
    # scatter 0.046mag
    V = G + 0.01760 + 0.006860*(Bp-Rp) + 0.1732*(Bp-Rp)**2

    #Please note that proper motion of RA is in seconds of time per year, as
    #required by Magellan telescope control software, which means: divide the
    #proper motion in arcsec by 15, then divide by cosine of Declination.
    # pmDEC: arcsec/yr
    pmDC = pmdec / 1e3
    pmDECstr = f'{pmDC:+.3f}'

    # pmRA: seconds/year
    pmRA = (pmra/(1e3*15))/np.cos(np.deg2rad(dec))
    pmRAstr = f'{pmRA:+.4f}'

    # last 15 characters of the EDR3 source_id become the target identifier
    # for this list
    assert len(source_id)>16
    targetname = source_id[-15:]

    if whichgaia in ['gaiadr2', 'gaiaedr3', 'gaiadr3']:
        equinox = '2000.0'
    else:
        raise NotImplementedError

    rowstr = (
        f'{targetname} {rahmsstr} {decdmsstr} {equinox} '+
        f'pmra={pmRAstr} pmdec={pmDECstr} vmag={V:.1f}'
    )

    return rowstr

def main():

    df = pd.read_csv("20210805_gmag_lt10_kinematicandpositionselected_backup_stars.csv")
    df = df.sort_values(by='ra')

    lines = []
    for ix,r in df.iterrows():
        s = str(r['source_id'])
        print(ix,s)
        r = get_log_row(s)
        lines.append(r+'\n')

    outpath = '20210806_hires_backup_targetlist.txt'
    with open(outpath, 'w') as f:
        f.writelines(lines)

if __name__ == "__main__":
    main()
