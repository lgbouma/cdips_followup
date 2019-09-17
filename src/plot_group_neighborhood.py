"""
You have a group. (Given by some large number of DR2 source_ids).

And a target star. (Usually one that shows transits, or generally one for which
you want to know if it is a group member).

Given these things, collect the neighbor stars. And make diagnostic plots that
let you know if it is or is not a convincing group member.
"""

###########
# imports #
###########
import os
import numpy as np, matplotlib.pyplot as plt, pandas as pd

from astropy.io.votable import from_table, writeto, parse
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

from astroquery.gaia import Gaia

from astrobase.services.gaia import objectid_search

datadir = '../data/'

homedir='/home/luke/'
credentials_file = os.path.join(homedir, '.gaia_credentials')
if not os.path.exists(credentials_file):
    raise AssertionError(
        'need gaia dr2 credentials file at {}'.format(credentials_file)
    )

gaiadir = os.path.join(homedir, '.gaia_cache')
if not os.path.exists(gaiadir):
    os.mkdir(gaiadir)


#############
# functions #
#############
def make_votable_given_source_ids(source_ids, outpath=None):

    t = Table()
    t['source_id'] = source_ids

    votable = from_table(t)

    writeto(votable, outpath)
    print('made {}'.format(outpath))

    return outpath


def given_votable_get_df(votablepath, assert_equal='source_id'):

    vot = parse(votablepath)
    tab = vot.get_first_table().to_table()
    df = tab.to_pandas()

    if isinstance(assert_equal, str):
        np.testing.assert_array_equal(tab[assert_equal], df[assert_equal])

    return df

def given_source_ids_get_gaia_data(source_ids, groupname, overwrite=True):
    """
    Args:

        source_ids (np.ndarray) of np.int64 Gaia DR2 source_ids

        groupname (str)

        overwrite: if True, and finds that this crossmatch has already run,
        deletes previous cached output and reruns anyway.

    Returns:

        dataframe with Gaia DR2 crossmatch info.
    """

    if type(source_ids) != np.ndarray:
        raise TypeError(
            'source_ids must be np.ndarray of np.int64 Gaia DR2 source_ids'
        )
    if type(source_ids[0]) != np.int64:
        raise TypeError(
            'source_ids must be np.ndarray of np.int64 Gaia DR2 source_ids'
        )

    xmltouploadpath = os.path.join(gaiadir,'toupload_{}.xml'.format(groupname))
    dlpath = os.path.join(gaiadir,'group{}_matches.xml.gz'.format(groupname))

    if overwrite:
        if os.path.exists(xmltouploadpath):
            os.remove(xmltouploadpath)

    if not os.path.exists(xmltouploadpath):
        make_votable_given_source_ids(source_ids, outpath=xmltouploadpath)

    if os.path.exists(dlpath) and overwrite:
        os.remove(dlpath)

    if not os.path.exists(dlpath):

        Gaia.login(credentials_file=credentials_file)

        jobstr = (
        '''
        SELECT *
        FROM tap_upload.foobar as u, gaiadr2.gaia_source AS g
        WHERE u.source_id=g.source_id
        '''
        )
        query = jobstr

        # might do async if this times out. but it doesn't.
        j = Gaia.launch_job(query=query,
                            upload_resource=xmltouploadpath,
                            upload_table_name="foobar", verbose=True,
                            dump_to_file=True, output_file=dlpath)

        Gaia.logout()

    df = given_votable_get_df(dlpath, assert_equal='source_id')

    if len(df) != len(source_ids):
        errmsg = (
            'ERROR! got {} matches vs {} source id queries'.
            format(len(df), len(source_ids))
        )
        print(errmsg)
        import IPython; IPython.embed()
        raise AssertionError(errmsg)

    return df


def query_neighborhood(bounds, groupname, n_max=2000, overwrite=True,
                       is_cg18_group=True, is_kc19_group=False):
    """
    Given the bounds in position and parallx corresponding to some group (e.g.,
    from Cantat-Gaudin+2018 or Kounkel & Covey 2019), get the DR2 stars from
    the group's neighborhood.

    The bounds are lower and upper in ra, dec, parallax, and there is a
    limiting G magnitude. A maximum number of stars, `n_max`, are selected from
    within these bounds.

    Args:
        bounds (dict): with parallax, ra, dec bounds.

        groupname (str)

        n_max (int): maximum number of stars in the neighborhood to acquire.

        is_kc19_group (bool): if the group is from Kounkel & Covey 2019,
        slightly different query is required.

    Returns:
        dataframe of DR2 stars within the bounds given. This is useful for
        querying stars that are in the neighborhood of some group.
    """

    if is_cg18_group:
        g_mag_limit=18
    elif is_kc19_group:
        g_mag_limit=18

    dlpath = os.path.join(
        gaiadir,'nbhd_group{}_matches.xml.gz'.format(groupname)
    )

    if os.path.exists(dlpath) and overwrite:
        os.remove(dlpath)

    if not os.path.exists(dlpath):

        Gaia.login(credentials_file=credentials_file)

        jobstr = (
        """
        select top {n_max:d}
            g.source_id, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
            g.phot_g_mean_mag, g.parallax, g.ra, g.dec, g.pmra, g.pmdec,
            g.radial_velocity, g.radial_velocity_error
        from gaiadr2.gaia_source as g
        where
            g.parallax > {parallax_lower:.2f}
        and
            g.parallax < {parallax_upper:.2f}
        and
            g.dec < {dec_upper:.2f}
        and
            g.dec > {dec_lower:.2f}
        and
            g.ra > {ra_lower:.2f}
        and
            g.ra < {ra_upper:.2f}
        and
            g.phot_g_mean_mag < {g_mag_limit:d}
        order by
            random_index
        """
        )
        query = jobstr.format(
            n_max=n_max,
            parallax_lower=bounds['parallax_lower'],
            parallax_upper=bounds['parallax_upper'],
            ra_lower=bounds['ra_lower'],
            ra_upper=bounds['ra_upper'],
            dec_lower=bounds['dec_lower'],
            dec_upper=bounds['dec_upper'],
            g_mag_limit=g_mag_limit
        )

        if is_kc19_group:
            # Kounkel & Covey impose some extra quality cuts.
            jobstr = (
            """
            select top {n_max:d}
                g.source_id, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
                g.phot_g_mean_mag, g.parallax, g.ra, g.dec, g.pmra, g.pmdec,
                g.radial_velocity, g.radial_velocity_error
            from gaiadr2.gaia_source as g
            where
                g.parallax > {parallax_lower:.2f}
            and
                g.parallax < {parallax_upper:.2f}
            and
                g.dec < {dec_upper:.2f}
            and
                g.dec > {dec_lower:.2f}
            and
                g.ra > {ra_lower:.2f}
            and
                g.ra < {ra_upper:.2f}
            and
                g.phot_g_mean_mag < {g_mag_limit:d}
            and
                parallax > 1
            and
                parallax_error < 0.1
            and
                1.0857/phot_g_mean_flux_over_error < 0.03
            and
                astrometric_sigma5d_max < 0.3
            and
                visibility_periods_used > 8
            and (
                    (astrometric_excess_noise < 1)
                    or
                    (astrometric_excess_noise > 1 and astrometric_excess_noise_sig < 2)
            )
            order by
                random_index
            """
            )
            query = jobstr.format(
                n_max=n_max,
                parallax_lower=bounds['parallax_lower'],
                parallax_upper=bounds['parallax_upper'],
                ra_lower=bounds['ra_lower'],
                ra_upper=bounds['ra_upper'],
                dec_lower=bounds['dec_lower'],
                dec_upper=bounds['dec_upper'],
                g_mag_limit=g_mag_limit
            )

        # async jobs can avoid timeout
        j = Gaia.launch_job_async(query=query, verbose=True, dump_to_file=True,
                                  output_file=dlpath)
        #j = Gaia.launch_job(query=query, verbose=True, dump_to_file=True,
        #                    output_file=dlpath)

        Gaia.logout()

    df = given_votable_get_df(dlpath, assert_equal='source_id')

    return df

def plot_group_neighborhood(
    targetname, groupname, group_df_dr2, target_df, nbhd_df,
    cutoff_probability,
    extra_overplot=0,
    pmdec_min=None,
    pmdec_max=None,
    pmra_min=None,
    pmra_max=None
):

    fig, axs = plt.subplots(2, 3, figsize=(18,12))

    ######################
    # parallax vs pm dec #
    ######################
    ax0 = axs[0,0]

    ax0.scatter(
        nbhd_df['pmdec'], nbhd_df['parallax'], c='gray', alpha=0.9, zorder=2,
        s=15, rasterized=True, linewidths=0, label='nbhd stars'
    )
    if extra_overplot:
        ax0.scatter(
            group_df_dr2['pmdec'], group_df_dr2['parallax'],
            c='k', alpha=0.9, zorder=3, s=15,
            rasterized=True, linewidths=0,
            label='CG18 P>{}'.format(cutoff_probability)
        )

    ax0.plot(
        target_df['pmdec'],
        target_df['parallax'],
        alpha=1, mew=0.5, zorder=8, label=targetname,
        markerfacecolor='yellow', markersize=16, marker='*', color='black',
        lw=0
    )

    ax0.legend(loc='best')
    ax0.set_xlabel(r'pmDEC, $\mu_{{\delta}}$ [mas/yr]', fontsize='xx-large')
    ax0.set_ylabel('star parallax [mas]', fontsize='xx-large')

    ax0.set_xlim([pmdec_min, pmdec_max])

    ######################
    # parallax vs pm ra #
    ######################
    ax1 = axs[0,1]

    ax1.scatter(
        nbhd_df['pmra'], nbhd_df['parallax'], c='gray', alpha=0.9, zorder=2, s=15,
        rasterized=True, linewidths=0, label='nbhd stars'
    )
    if extra_overplot:
        ax1.scatter(
            group_df_dr2['pmra'], group_df_dr2['parallax'], c='k', alpha=0.9,
            zorder=3, s=15, rasterized=True, linewidths=0,
            label='CG18 P>{}'.format(cutoff_probability)
        )
    ax1.plot(
        target_df['pmra'], target_df['parallax'], alpha=1, mew=0.5, zorder=8,
        label=targetname, markerfacecolor='yellow', markersize=16, marker='*',
        color='black', lw=0
    )

    ax1.legend(loc='best')
    ax1.set_xlabel(r'pmRA, $\mu_{{\alpha}} \cos\delta$ [mas/yr]',
                   fontsize='xx-large')
    ax1.set_ylabel('star parallax [mas]', fontsize='xx-large')

    ax1.set_xlim([pmra_min, pmra_max])

    ##################
    # proper motions #
    ##################
    ax2 = axs[0,2]

    ax2.scatter(
        nbhd_df['pmra'], nbhd_df['pmdec'], c='gray', alpha=0.9, zorder=2, s=15,
        rasterized=True, linewidths=0, label='nbhd stars'
    )
    if extra_overplot:
        ax2.scatter(
            group_df_dr2['pmra'], group_df_dr2['pmdec'], c='k', alpha=0.9, zorder=3, s=15,
            rasterized=True, linewidths=0,
            label='CG18 P>{}'.format(cutoff_probability)
        )
    ax2.plot(
        target_df['pmra'], target_df['pmdec'], alpha=1, mew=0.5, zorder=8,
        label=targetname, markerfacecolor='yellow', markersize=16, marker='*',
        color='black', lw=0
    )

    ax2.legend(loc='best')

    ax2.set_xlabel(r'pmRA, $\mu_{{\alpha}} \cos\delta$ [mas/yr]',
                   fontsize='xx-large')
    ax2.set_ylabel(r'pmDEC, $\mu_{{\delta}}$ [mas/yr]',
                   fontsize='xx-large')

    ax2.set_xlim([pmra_min, pmra_max])
    ax2.set_ylim([pmdec_min, pmdec_max])

    ##############
    # HR diagram # 
    ##############
    ax3 = axs[1,0]

    ax3.scatter(
        nbhd_df['phot_bp_mean_mag']-nbhd_df['phot_rp_mean_mag'], nbhd_df['phot_g_mean_mag'],
        c='gray', alpha=1., zorder=2, s=15, rasterized=True, linewidths=0,
        label='nbhd stars'
    )
    ax3.scatter(
        group_df_dr2['phot_bp_mean_mag']-group_df_dr2['phot_rp_mean_mag'],
        group_df_dr2['phot_g_mean_mag'],
        c='k', alpha=1., zorder=3, s=15, rasterized=True, linewidths=0,
        label='CG18 P>{}'.format(cutoff_probability)
    )
    ax3.plot(
        target_df['phot_bp_mean_mag']-target_df['phot_rp_mean_mag'],
        target_df['phot_g_mean_mag'],
        alpha=1, mew=0.5, zorder=8, label=targetname, markerfacecolor='yellow',
        markersize=16, marker='*', color='black', lw=0
    )

    ax3.legend(loc='best')
    ax3.set_ylabel('G', fontsize='xx-large')
    ax3.set_xlabel('Bp - Rp', fontsize='xx-large')

    ylim = ax3.get_ylim()
    ax3.set_ylim((max(ylim),min(ylim)))
    ax3.set_xlim((-0.5, 3.0))

    #############
    # positions #
    #############
    ax4 = axs[1,1]

    ax4.scatter(
        nbhd_df['ra'], nbhd_df['dec'], c='gray', alpha=1., zorder=2, s=15,
        rasterized=True, linewidths=0, label='nbhd stars'
    )
    ax4.scatter(
        group_df_dr2['ra'], group_df_dr2['dec'], c='k', alpha=1., zorder=3,
        s=15, rasterized=True, linewidths=0,
        label='CG18 P>{}'.format(cutoff_probability)
    )
    ax4.plot(
        target_df['ra'], target_df['dec'], alpha=1, mew=0.5, zorder=8,
        label=targetname, markerfacecolor='yellow', markersize=16, marker='*',
        color='black', lw=0
    )

    ax4.legend(loc='best')
    ax4.set_xlabel(r'RA, $\alpha$ [deg]', fontsize='xx-large')
    ax4.set_ylabel('Dec, $\delta$ [deg]', fontsize='xx-large')

    ######################
    # parallax vs DR2 rv #
    ######################
    ax5 = axs[1,2]

    ax5.scatter(
        nbhd_df['radial_velocity'], nbhd_df['parallax'], c='gray', alpha=0.9,
        zorder=2, s=15, rasterized=True, linewidths=0, label='nbhd stars'
    )
    ax5.plot(
        target_df['radial_velocity'], target_df['parallax'], alpha=1, mew=0.5,
        zorder=8, label=targetname, markerfacecolor='yellow', markersize=16,
        marker='*', color='black', lw=0
    )

    ax5.legend(loc='best')
    ax5.set_xlabel('DR2 RV (if available) [km/s]', fontsize='xx-large')
    ax5.set_ylabel('star parallax [mas]', fontsize='xx-large')


    ########
    # save #
    ########
    for ax in [ax0,ax1,ax2,ax3,ax4,ax5]:
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='xx-large')
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='xx-large')

    fig.tight_layout()

    outpath = '../results/{}_neighborhood.png'.format(groupname)
    if extra_overplot:
        outpath = '../results/{}_neighborhood_extra.png'.format(groupname)
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print('made {}'.format(outpath))


def main(extra_overplot=0):

    overwrite = 1
    groupname = 'NGC2354'
    targetname = '554-051'
    target_id = 5617126180115568256
    cutoff_probability = 0.1

    ##########################################
    # some data for a group that have a source_id list
    cg18_path = os.path.join(
        datadir, 'ngc_2354_CG18_subset.vot.gz')
    group_df = given_votable_get_df(cg18_path, assert_equal='Source')
    group_df = group_df[group_df['PMemb'] > cutoff_probability]
    group_source_ids = np.array(group_df['Source']).astype(np.int64)
    # todo: generalize above -- generically, it just needs some source id list.
    ##########################################

    group_df_dr2 = given_source_ids_get_gaia_data(group_source_ids, groupname,
                                                  overwrite=overwrite)

    target_d = objectid_search(target_id)
    target_df = pd.read_csv(target_d['result'])
    assert len(target_df) == 1

    # now acquire the mean properties of the group, and query the neighborhood
    # based on those properties. the number of neighbor stars to randomly
    # select is min(5* the number of group members, 5000).
    bounds = {}
    params = ['parallax', 'ra', 'dec']
    for param in params:

        bounds[param+'_upper'] = (
            group_df_dr2[param].mean() + 5*group_df_dr2[param].std()
        )

        bounds[param+'_lower'] = (
            group_df_dr2[param].mean() - 5*group_df_dr2[param].std()
        )

    assert bounds['ra_upper'] < 360
    assert bounds['ra_lower'] > 0

    n_max = min((50*len(group_df_dr2), 10000))
    nbhd_df = query_neighborhood(bounds, groupname, n_max=n_max,
                                 overwrite=overwrite, is_cg18_group=True,
                                 is_kc19_group=False)

    # ensure no overlap between the group members and the neighborhood sample.
    common = group_df_dr2.merge(nbhd_df, on='source_id', how='inner')
    snbhd_df = nbhd_df[~nbhd_df.source_id.isin(common.source_id)]

    plot_group_neighborhood(targetname, groupname, group_df_dr2, target_df,
                            nbhd_df, cutoff_probability,
                            extra_overplot=extra_overplot)



if __name__ == "__main__":
    main(extra_overplot=0)
    main(extra_overplot=1)
