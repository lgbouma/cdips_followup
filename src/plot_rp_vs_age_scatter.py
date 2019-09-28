'''
DESCRIPTION
----------
Download exoplanet archive and make scatter plot of Rp vs age.

* R_p vs age, density vs age
    -> are hot jupiters typically around younger stars?
    -> can we see photoevaporation in time? [less dense planets lose
       atmospheres, get smaller and more dense]
    -> do we see the radius gap move in time?

USAGE
----------
select desired plots from bools below. then:

$ python exoarchive_age_plots.py

'''

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u

import os

from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

def arr(x):
    return np.array(x)


if __name__=='__main__':
    show_median_error=1

    # columns described at
    # https://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html
    ea_tab = NasaExoplanetArchive.get_confirmed_planets_table(all_columns=True,
                                                              show_progress=True)

    # radius in jupiter (earth) radii, period, number of planets in system,
    # galacitc latitude.
    wellmeasuredparams = ['pl_radj', 'pl_rade', 'pl_orbper', 'pl_pnum', 'st_glat']
    # janky measurements where errors matter: age in Gyr, density(g/cc), eccen,
    # inclination
    jankyparams = ['pl_dens', 'pl_orbeccen', 'pl_orbincl']

    # get systems with finite ages (has a value, and +/- error bar)
    has_age_value = ~ea_tab['st_age'].mask
    has_age_errs  = (~ea_tab['st_ageerr1'].mask) & (~ea_tab['st_ageerr2'].mask)
    has_rp_value = ~ea_tab['pl_rade'].mask
    has_rp_errs  = (~ea_tab['pl_radeerr1'].mask) & (~ea_tab['pl_radeerr2'].mask)

    transits = (ea_tab['pl_tranflag']==1)

    sel = (
        has_age_value & has_age_errs & has_rp_value & has_rp_errs & transits
    )

    t = ea_tab[sel]

    #
    # read params
    #
    age = t['st_age']
    age_perr = t['st_ageerr1']
    age_merr = np.abs(t['st_ageerr2'])
    age_errs = np.array([age_perr, age_merr]).reshape(2, len(age))

    rp = t['pl_rade']
    rp_perr = t['pl_radeerr1']
    rp_merr = t['pl_radeerr2']
    rp_errs = np.array([rp_perr, rp_merr]).reshape(2, len(age))


    # plot age vs all "good cols". (age is on y axis b/c it has the error bars, and
    # I at least skimmed the footnotes of Hogg 2010)
    f,ax = plt.subplots(figsize=(4,3))

    # squares
    label = (
            'NASA exoarchive (median $\sigma_{{\mathrm{{age}}}}$ = {:.1f} Gyr)'.
            format(np.median(age_errs))
    )
    ax.scatter(age, rp,
               color='gray', s=2, zorder=1, marker='o', linewidth=0,
               label=label, alpha=0.6)

    # targets
    tdf = pd.read_csv(
        '/home/luke/Dropbox/proj/cdips_followup/results/20190926_2020A_targets_age_rp.csv'
    )
    target_age = 10**(np.array(tdf['age']))/(1e9)
    target_rp = np.array(tdf['rplanet'])
    target_rp_unc = np.array(tdf['rplanet_unc'])

    ax.scatter(
        target_age, target_rp, color='black', s=20,
        zorder=3, marker='*', linewidth=0,
        label='Proposed targets (median $\sigma_{{\mathrm{{age}}}}$ < 0.1 Gyr)'
    )

    ax.errorbar(target_age, target_rp, yerr=target_rp_unc,
                elinewidth=0.3, ecolor='k', capsize=0, capthick=0, linewidth=0,
                fmt='*', ms=0, zorder=2, color='black', alpha=0.5)

    ## error bars
    #ax.errorbar(age, rp, yerr=age_errs, xerr=rp_errs,
    #            elinewidth=0.3, ecolor='k', capsize=0, capthick=0, linewidth=0,
    #            fmt='s', ms=0, zorder=2, alpha=0.05)

    #if show_median_error:
    #    ax.text(0.03, 0.03, txt, ha='left', va='bottom', fontsize='small',
    #            zorder=3, transform=ax.transAxes, color='C0')

    ax.legend(loc='best', borderpad=0.1, handletextpad=0.1, fontsize='small',
              framealpha=1)

    ax.set_xlabel('Age [Gyr]')
    ax.set_ylabel('Planet radius [R$_\oplus$]')

    #ax.tick_params(top=True, bottom=True, left=True, right=True)
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax.set_xscale('log')
    ax.set_yscale('log')

    f.tight_layout()
    f.savefig('../results/rp_vs_age_scatter.png', bbox_inches='tight', dpi=300)
    f.savefig('../results/rp_vs_age_scatter.pdf', bbox_inches='tight')
