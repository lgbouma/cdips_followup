'''
DESCRIPTION
----------
Make scatter plot of Rp vs age using exoplanet archive and candidates.csv.
(Note: no uncertainties, because the Rp uncertainties aren't in the
candidates.csv database).

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

def plot_rp_vs_age_scatter(active_targets=0, split_toi_ctoi=0):

    #
    # columns described at
    # https://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html
    #
    ea_tab = NasaExoplanetArchive.get_confirmed_planets_table(
        all_columns=True, show_progress=True
    )

    #
    # get systems with finite ages (has a value, and +/- error bar)
    #
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

    #
    # plot age vs rp. (age is on y axis b/c it has the error bars, and I at
    # least skimmed the footnotes of Hogg 2010)
    #
    f,ax = plt.subplots(figsize=(4,3))

    label = (
        'NASA Exoplanet Archive, '+
        r"$\langle \sigma_{{\mathrm{{age}}}} \rangle$ = "+
        "{:.1f} Gyr".format(np.median(age_errs))
    )
    ax.scatter(age, rp,
               color='gray', s=3, zorder=1, marker='o', linewidth=0,
               label=label, alpha=1)

    #
    # targets
    #
    cdf = pd.read_csv(
        os.path.join(
            os.path.expanduser("~"),
            'Dropbox/proj/cdips_followup/data/candidate_database/candidates.csv'
        ),
        sep='|'
    )

    sel = (
        ~cdf.isretired
        &
        (cdf.current_priority <= 1)
        &
        ~pd.isnull(cdf.rp)
        &
        cdf.iscdipstarget
    )

    sdf = cdf[sel]

    target_age = np.array(sdf.age)
    target_rp = np.array(sdf.rp)
    target_rp_unc = np.array(sdf.rp_unc)

    temp_ages = []
    for a in target_age:

        if a == '--':
            temp_ages.append('--')
            continue

        temp_ages.append(
            np.mean(np.array(a.split(',')).astype(float))
        )

    target_age = np.array(temp_ages)

    #
    # fix non-assigned ages.
    # 1) "PMS" star with no age -> 500 Myr upper bound.
    # 2) Vela OB2 subgroups get ages according to Cantat-Gaudin2019, Figure 6.
    #    (Use the "name" column and match the "cg19velaOB2_pop[N]" pattern).
    #
    is_pms = (sdf.reference == 'Zari_2018_PMS') & (sdf.age == '--')
    target_age[is_pms] = np.log10(5e8)

    vela_ob2_age_dict = {
        '1': np.log10(5e7),
        '2': np.log10(4e7),
        '3': np.log10(4e7),
        '4': np.log10(3e7),
        '5': np.log10(3e7),
        '6': np.log10(2.5e7),
        '7': np.log10(1.5e7)
    }

    popn_inds = np.array(
        sdf.name.str.extract(pat='cg19velaOB2_pop(\d)')
    ).flatten()

    cg19_ages = []
    for k in popn_inds:
        if pd.isnull(k):
            cg19_ages.append(np.nan)
        else:
            cg19_ages.append(vela_ob2_age_dict[k[0]])

    target_age[~pd.isnull(cg19_ages)] = np.array(cg19_ages)[~pd.isnull(cg19_ages)]

    target_age = target_age.astype(float)

    target_age = 10**(np.array(target_age))/(1e9)

    if active_targets:

        label = (
            r'Active targets, $\langle \sigma_{{\mathrm{{age}}}} \rangle$'
            '< 0.1 Gyr'
        )

        istoi = ~(sdf.toi == '--')

        if split_toi_ctoi:
            ax.scatter(
                target_age[istoi], target_rp[istoi], color='black', s=10, zorder=3,
                marker='s', linewidth=0
            )
            ax.scatter(
                target_age[~istoi], target_rp[~istoi], color='black', s=25,
                zorder=3, marker='*', linewidth=0, label=label
            )
        else:
            ax.scatter(
                target_age, target_rp, color='black', s=25,
                zorder=3, marker='*', linewidth=0, label=label
            )

        target_rp_rel_unc = target_rp_unc / target_rp

        ##########
        # HACK: force to 30% relative uncertainty for TOI-520, b/c fits did not
        # converge
        sel = target_rp_rel_unc > 0.5
        target_rp_unc[sel] = target_rp[sel] * 0.3
        ##########

        ax.errorbar(target_age, target_rp, yerr=target_rp_unc, elinewidth=0.3,
                    ecolor='k', capsize=0, capthick=0, linewidth=0, fmt='*',
                    ms=0, zorder=2, color='black', alpha=0.5)

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
    if split_toi_ctoi:
        savstr += '_split_ctoi'

    outpath = '../results/rp_vs_age_scatter{}.png'.format(savstr)
    f.savefig(outpath, bbox_inches='tight', dpi=450)
    f.savefig(outpath.replace('.png','.pdf'), bbox_inches='tight')
    print('made {}'.format(outpath))

if __name__=='__main__':

    plot_rp_vs_age_scatter(active_targets=0)
    plot_rp_vs_age_scatter(active_targets=1, split_toi_ctoi=0)
    plot_rp_vs_age_scatter(active_targets=1, split_toi_ctoi=1)
