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
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter

import pandas as pd, numpy as np
import os
from copy import deepcopy

from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u

from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

from cdips.utils import today_YYYYMMDD
from numpy import array as arr
from aesthetic.plot import set_style

def plot_rp_vs_age_scatter(active_targets=0, split_toi_ctoi=0, hjs_only=0,
                           ismanualsubset=0, isvalidated=0, ispublictalk=0):

    #
    # columns described at
    # https://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html
    #
    ea_tab = NasaExoplanetArchive.query_criteria(
        table="exoplanets", select="*", cache=True
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

    if hjs_only:
        is_hj = (t['pl_orbper'] < 10*u.day) & (t['pl_rade'] > 7)
        t = t[is_hj]

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
    set_style()
    f,ax = plt.subplots(figsize=(4,3))

    label = (
        'Exoplanet Archive, '+
        r"$\langle \sigma_{{\mathrm{{age}}}} \rangle$ = "+
        "{:.1f} Gyr".format(np.median(age_errs))
    )
    color = 'C1'
    if hjs_only:
        label = (
            'Known hot Jupiters'
        )
        rp /= 11.2089 # used jupiter radii
        color = 'C1'
    if ispublictalk:
        label = 'Known'
        color = 'lightgray'

    ax.scatter(age, rp,
               color=color, s=3, zorder=1, marker='o', linewidth=0,
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

    if isvalidated:
        valsourceids = [
            "5251470948229949568", # TOI837
        ]
        val_df = pd.DataFrame(valsourceids, columns=['source_id'])
        cdf.source_id = cdf.source_id.astype(str)
        mdf = val_df.merge(cdf, on='source_id')
        vdf = deepcopy(mdf)
        assert len(vdf) == len(valsourceids)

    if ismanualsubset:
        manualsourceids = [
            # "5523449717870971776", # probable EB
            "5489726768531119616",
            "5510676828723793920",
            "1835201042675810688",
            "6598814657249555328",
            "5557593814516968960",
            "5952590785523816960",
            "5514577414951631488",
            "4844691297067063424",
            "5525188767305211904",
            "5838450865699668736",
            "5765748511163751936",
            # "6113920619134019456", # Rizzuto+20 thyme2
            "6042883578050870912",
            # "5974331982990013696", # suggestive, not good enough for PC
            "5519619186857962112",
            "5838183443852841216",
            "5239758155778687360",
            "5254512781523942912"
        ]
        source_df = pd.DataFrame(manualsourceids, columns=['source_id'])
        cdf.source_id = cdf.source_id.astype(str)
        mdf = source_df.merge(cdf, on='source_id')
        sdf = deepcopy(mdf)
        assert len(sdf) == len(manualsourceids)

    else:
        # standard selection
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
        if hjs_only:
            label = (
                'Possible hot Jupiters'
            )

        istoi = ~(sdf.toi == '--')

        if split_toi_ctoi:
            ax.scatter(
                target_age[istoi], target_rp[istoi], color='C0', s=10, zorder=3,
                marker='s', linewidth=0
            )
            ax.scatter(
                target_age[~istoi], target_rp[~istoi], color='C0', s=25,
                zorder=3, marker='*', linewidth=0, label=label
            )
        elif hjs_only:
            hj = (target_rp > 7) & (target_rp < 29)
            ax.scatter(
                target_age[hj], target_rp[hj]/11.2089, color='C0', s=35,
                zorder=3, marker='*', linewidth=0, label=label
            )
        else:
            if isvalidated:
                if len(vdf) > 1:
                    raise NotImplementedError('might want to just assign ages')
                val_age = 4e7/1e9
                val_rp = np.array(vdf.rp)

                ax.plot(
                    val_age, val_rp, mew=0.5, markerfacecolor='yellow',
                    markersize=12, marker='*', color='k', lw=0,
                    label='Validated', zorder=4
                )

                ax.scatter(
                    target_age, target_rp, s=40, color='C0',
                    zorder=3, marker='*', linewidth=0, label="Potential"
                )

            else:
                ax.scatter(
                    target_age, target_rp, color='C0', s=25,
                    zorder=3, marker='*', linewidth=0, label=label
                )

        sel = (~pd.isnull(target_age)) & (~pd.isnull(target_rp))
        print('{} targets w/ finite ages and radii'.
              format(len(target_age[sel])))
        print('{} targets w/ finite ages and radii that are TOIs'.
              format(len(target_age[sel & istoi])))
        print('{} targets w/ finite ages and radii that are not TOIs'.
              format(len(target_age[sel & ~istoi])))

        target_rp_rel_unc = target_rp_unc / target_rp

        ##########
        if hjs_only:
            # HACK: force to 10% relative uncertainty (eep)
            sel = target_rp_rel_unc > 0.1
            target_rp_unc[sel] = target_rp[sel] * 0.1 * (
                np.random.uniform(0.8,1.2,
                                  size=len(target_rp[sel]))
            )
        else:
            # HACK: for the one very small object...
            sel = target_rp_rel_unc > 0.33
            target_rp_unc[sel] = target_rp[sel] * 0.33
        ##########

        if hjs_only:
            pass
            #ax.errorbar(target_age[hj], target_rp[hj], yerr=target_rp_unc[hj],
            #            elinewidth=0.3, ecolor='k', capsize=0, capthick=0,
            #            linewidth=0, fmt='*', ms=0, zorder=2, color='black',
            #            alpha=0.5)
        else:
            if not ispublictalk:
                ax.errorbar(target_age, target_rp, yerr=target_rp_unc,
                            elinewidth=0.3, ecolor='k', capsize=0, capthick=0,
                            linewidth=0, fmt='*', ms=0, zorder=2, color='black',
                            alpha=0.5)

    # flip default legend order
    handles, labels = ax.get_legend_handles_labels()
    if not hjs_only:
        leg = ax.legend(handles[::-1], labels[::-1], loc='lower left',
                        borderpad=0.3, handletextpad=0.3, fontsize='small',
                        framealpha=1)
    else:
        leg = ax.legend(handles[::-1], labels[::-1], loc='upper right',
                        borderpad=0.3, handletextpad=0.3, fontsize='small',
                        framealpha=0)

    leg.get_frame().set_linewidth(0.5)

    ax.set_xlabel('Planet age [billion years]')
    ax.set_ylabel('Planet size [Earth radii]')
    if hjs_only:
        ax.set_ylabel('Planet size [Jupiter radii]')

    ax.set_xlim([6e-3, 17])

    if not hjs_only:
        ax.set_ylim([0.13, 85])
        ax.set_yscale('log')
    if hjs_only:
        ax.set_ylim([0.1, 3])
        ax.set_yscale('linear')

    ax.set_xscale('log')

    if hjs_only:
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.2g}'))
    else:
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.2g}'))


    f.tight_layout()
    savstr = '_no_overplot' if not active_targets else '_active_targets'
    if split_toi_ctoi:
        savstr += '_split_ctoi'
    if hjs_only:
        savstr += '_onlyhjs'
    if ispublictalk:
        savstr += '_ispublic'
    if isvalidated:
        savstr += '_hasvalidated'

    outpath = (
        '../results/rp_vs_age_scatter/rp_vs_age_scatter_{}{}.png'.
        format(today_YYYYMMDD(), savstr)
    )
    f.savefig(outpath, bbox_inches='tight', dpi=450)
    f.savefig(outpath.replace('.png','.pdf'), bbox_inches='tight')
    print('made {}'.format(outpath))


if __name__=='__main__':

    active_targets = 1  # show current follow-up targets.
    ismanualsubset = 1  # if false, will show all current follow-up targets.
    isvalidated = 1     # whether to include CDIPS paper as their own thing.
    split_toi_ctoi = 0  # whether to split based on TOI/CTOI
    hjs_only = 0        # changes units to jupiter radii
    ispublictalk = 1    # changes labels and units

    plot_rp_vs_age_scatter(active_targets=active_targets,
                           split_toi_ctoi=split_toi_ctoi,
                           ismanualsubset=ismanualsubset,
                           hjs_only=hjs_only,
                           ispublictalk=ispublictalk,
                           isvalidated=isvalidated)
