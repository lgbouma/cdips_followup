"""
%\begin{targettable}{}
%\objid{}           % specify a 3-digit number for each target.
%\object{}          % 20 characters maximum
%\ra{}              % e.g., xx:xx:xx.x HH:MM:SS.S
%\dec{}             % e.g., +-xx:xx:xx.x
%\epoch{}           % e.g., 1950.3
%\magnitude{}
%\filter{}
%\exptime{}         % in seconds PER EXPOSURE
%\nexposures{}      % Number of exposures
%\moondays{}        % Days from new moon, use a number 0-14
%\skycond{}         % "spec" or "phot"
%\seeing{}          % max allowable PSF FWHM (arcsecs)
%\obscomment{}      % 20 characters maximum - REQUIRED COMMAND
%  - repeat target entry parameters as needed to complete Table -
%\end{targettable}
"""

from glob import glob
import os
import pandas as pd, numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c
from pfs_calculator import given_Vmag_get_PFS_exptime_minutes

def main():

    df = pd.read_csv('../results/20190912_2020A_targets_obs_info.csv', sep='|')

    # TOI450, and 4 objects that are just too faint (gaia G > 15).
    # And no Teff>8000K objects
    badids = [
    4827527233363019776, # TOI 450
    2919143383943171200, # 4 objects that have G>15, too faint
    3340674976430098688,
    5561614350584396800,
    5618515825371166464,
    5519619186857962112, # No Teff > 8000K objects
    5525188767305211904,
    5290721997195236480,
    5557593814516968960,
    5579734916388215808,
    5325454783543157760,
    5516140233292943872 # PMS M dwarf (too faint in V)
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
