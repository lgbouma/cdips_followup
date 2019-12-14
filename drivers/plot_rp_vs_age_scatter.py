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

def plot_rp_vs_age_scatter(active_targets=0):

    # columns described at
    # https://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html
    ea_tab = NasaExoplanetArchive.get_confirmed_planets_table(
        all_columns=True, show_progress=True
    )

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
        'NASA Exoplanet Archive, '+
        r"$\langle \sigma_{{\mathrm{{age}}}} \rangle$ = "+
        "{:.1f} Gyr".format(np.median(age_errs))
    )
    ax.scatter(age, rp,
               color='gray', s=3, zorder=1, marker='o', linewidth=0,
               label=label, alpha=1)

    # targets
    tdf = pd.read_csv(
        os.path.join(
            os.path.expanduser("~"),
            'Dropbox/proj/cdips_followup/results/20190926_2020A_targets_age_rp.csv'
        )
    )

    badids = [
    4827527233363019776,
    2919143383943171200,
    3340674976430098688,
    5561614350584396800,
    5618515825371166464
    ]
    tdf = tdf[~tdf['source_id'].isin(badids)]

    target_age = 10**(np.array(tdf['age']))/(1e9)
    target_rp = np.array(tdf['rplanet'])
    target_rp_unc = np.array(tdf['rplanet_unc'])

    if active_targets:
        ax.scatter(
            target_age, target_rp, color='black', s=25,
            zorder=3, marker='*', linewidth=0,
            label=r'Active targets, $\langle \sigma_{{\mathrm{{age}}}} \rangle$ < 0.1 Gyr'
        )

        target_rp_rel_unc = target_rp_unc / target_rp

        ##########
        # HACK: force to 30% relative uncertainty for TOI-520, b/c fits did not
        # converge
        sel = target_rp_rel_unc > 0.5
        target_rp_unc[sel] = target_rp[sel] * 0.3
        ##########

        ax.errorbar(target_age, target_rp, yerr=target_rp_unc,
                    elinewidth=0.3, ecolor='k', capsize=0, capthick=0, linewidth=0,
                    fmt='*', ms=0, zorder=2, color='black', alpha=0.5)

    ## error bars
    #ax.errorbar(age, rp, yerr=age_errs, xerr=rp_errs,
    #            elinewidth=0.3, ecolor='k', capsize=0, capthick=0, linewidth=0,
    #            fmt='s', ms=0, zorder=2, alpha=0.05)

    # flip default legend order
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='lower left', borderpad=0.1,
              handletextpad=0.1, fontsize='small', framealpha=1)

    ax.set_xlabel('Age [Gyr]')
    ax.set_ylabel('Planet radius [R$_\oplus$]')

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax.set_ylim([0.13, 85])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(top=True, bottom=True, left=True, right=True, which='both')


    f.tight_layout()
    savstr = '_no_overplot' if not active_targets else '_active_targets'

    outpath = '../results/rp_vs_age_scatter{}.png'.format(savstr)
    f.savefig(outpath, bbox_inches='tight', dpi=450)
    f.savefig(outpath.replace('.png','.pdf'), bbox_inches='tight')
    print('made {}'.format(outpath))

if __name__=='__main__':
    plot_rp_vs_age_scatter(active_targets=0)
    plot_rp_vs_age_scatter(active_targets=1)
