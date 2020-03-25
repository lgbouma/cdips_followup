"""
DEPRECATED

I realized that a google sheet view and screenshot are 1000x better.
"""

from glob import glob
import os
import pandas as pd, numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c
from cdips_followup.exposure_calculators.pfs_calculator import given_Vmag_get_PFS_exptime_minutes

def main():

    # targes that were chosen from CDIPS_candidates_read-only after
    # applying relevant PRV filters (Teff<6800K, priority 1 or 0,
    # G<13.5, no SP2).
    df = pd.read_csv('../../results/followup_target_selection/20200325_2020B_target_ids.csv')

    # can list bad sourceids here.
    badids = [
    ]

    sdf = df[~df['source_id'].isin(badids)]

    objid = 1

    txtrows = []

    for ix, r in sdf.iterrows():

        c = SkyCoord(float(r['ra']), float(r['dec']), frame='icrs',
                     unit=(u.deg, u.deg))

        epoch = 2015.5

        _objid = "\\objid{{{objid:s}}}\n".format(
            objid=str(objid).zfill(3))

        _object = "\\object{{{_object:s}}}\n".format(
            _object=str(r['toi_or_ticid']).replace('TIC',''))

        ra = "\\ra{{{ra:s}}}\n".format(
            ra=c.ra.to_string(unit=u.hourangle, precision=1, sep=':'))

        dec = "\\dec{{{dec:s}}}\n".format(
            dec=c.dec.to_string(unit=u.deg, precision=1, sep=':', alwayssign=True))

        epoch = "\\epoch{{{epoch:.1f}}}\n".format(
            epoch=epoch)

        gmag = float(r['phot_g_mean_mag'])
        magnitude = "\\magnitude{{{magnitude:.1f}}}\n".format(
            magnitude=gmag
        )

        mode = '-'
        _filter = "\\filter{{{mode:s}}}\n".format(mode=mode.replace('_mode',''))

        exptime = given_Vmag_get_PFS_exptime_minutes(gmag)

        exptime = np.int(exptime)
        if exptime == 0:
            exptime = 1

        _exptime = "\\exptime{{{exptime:d}}}\n".format(
            exptime=exptime
        )

        nexposures = '-'

        _nexposures = "\\nexposures{{{nexposures:s}}}\n".format(
            nexposures=nexposures
        )

        if gmag < 10.5:
            moondays = 14
        elif gmag < 14:
            moondays = 11
        elif gmag < 16:
            moondays = 11
        moondays = "\\moondays{{{moondays:d}}}\n".format(
            moondays=moondays
        )

        skycond = "\\skycond{-}\n"
        seeing = "\\seeing{5}\n"
        obscomment = "\\obscomment{{{source_id:s}}}".format(
            source_id=str(r['source_id'])
        )

        # NOTE: deprecated
        txtrow = (
            _objid+_object+ra+dec+epoch+magnitude+_filter+
            _exptime+_nexposures+moondays+skycond+seeing+obscomment+'\n%\n'
        )

        txtrows.append(txtrow)

        objid += 1

    outlines = ''.join(txtrows)
    outpath = '../results/2020A_PFS_target_table.tex'
    with open(outpath, 'w') as f:
        f.writelines(outlines)
    print('wrote {}'.format(outpath))

if __name__ == "__main__":
    main()
