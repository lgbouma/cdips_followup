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
from cdips_followup.exposure_calculators.chiron_calculator import (
    get_chiron_reqd_exptime,
    get_chiron_snr_given_exptime
)

def main():

    df = pd.read_csv('../results/20190912_2020A_targets_obs_info.csv', sep='|')

    # TOI450, and 4 objects that are just too faint (gaia G > 15).
    badids = [
    4827527233363019776,
    2919143383943171200,
    3340674976430098688,
    5561614350584396800,
    5618515825371166464
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

        if gmag < 10.5:
            mode = 'slicer_mode'
        else:
            mode = 'fiber_mode'
        _filter = "\\filter{{{mode:s}}}\n".format(mode=mode.replace('_mode',''))

        exptime = get_chiron_reqd_exptime(gmag, mode,
                                          target_snr_per_resoln_element_all_exposures=50,
                                          n_exposures=1, verbose=False)

        # round up to the nearest 10
        exptime = np.int(np.round(exptime+10, decimals=-1))

        _exptime = "\\exptime{{{exptime:d}}}\n".format(
            exptime=exptime
        )
        if exptime > 1800:

            print('WARNING: GOT EXPTIME > 30 MIN, {:d}sec, for {}, Gmag {:.1f}'.
                  format(exptime, str(r['toi_or_ticid']), gmag))

            wouldgetsnr = get_chiron_snr_given_exptime(
                gmag, mode, exptime, verbose=False)

            willgetsnr = get_chiron_snr_given_exptime(
                gmag, mode, 1800, verbose=False)

            _exptime = "\\exptime{{{exptime:d}}}\n".format(
                exptime=1800
            )

            print('\tWould get SNR = {:.1f}. At 1800sec, will get SNR = {:.1f}'.
                  format(wouldgetsnr, willgetsnr))

        nexposures = 3
        #if 3600 > exptime > 1800:
        #    nexposures = 6
        #    _exptime = "\\exptime{{{exptime:d}}}\n".format(
        #        exptime=1800
        #    )

        #if 6100 > exptime > 3600:
        #    # 6100 sec is the faintest 2020A target. mainly a case-study to see
        #    # if this works.
        #    nexposures = 9
        #    _exptime = "\\exptime{{{exptime:d}}}\n".format(
        #        exptime=1800
        #    )

        _nexposures = "\\nexposures{{{nexposures:d}}}\n".format(
            nexposures=nexposures
        )

        if gmag < 10.5:
            moondays = 14
        elif gmag < 14:
            moondays = 11
        elif gmag < 16:
            moondays = 7
        moondays = "\\moondays{{{moondays:d}}}\n".format(
            moondays=moondays
        )

        skycond = "\\skycond{spec}\n"
        seeing = "\\seeing{5}\n"
        obscomment = "\\obscomment{{{source_id:s}}}".format(
            source_id=str(r['source_id'])
        )
        if nexposures == 6:
            obscomment = "\\obscomment{{2 exp per visit (x3)}}"
        if nexposures == 9:
            obscomment = "\\obscomment{{3 exp per visit (x3)}}"

        txtrow = (
            _objid+_object+ra+dec+epoch+magnitude+_filter+
            _exptime+_nexposures+moondays+skycond+seeing+obscomment+'\n%\n'
        )

        txtrows.append(txtrow)

        objid += 1

    outlines = ''.join(txtrows)
    outpath = '../results/2020A_NOAO_target_table.tex'
    with open(outpath, 'w') as f:
        f.writelines(outlines)
    print('wrote {}'.format(outpath))

if __name__ == "__main__":
    main()
