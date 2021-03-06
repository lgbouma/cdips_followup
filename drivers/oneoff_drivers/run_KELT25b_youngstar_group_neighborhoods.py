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
from astroquery.vizier import Vizier

from astrobase.services.gaia import objectid_search

from plot_group_neighborhood import (
    make_votable_given_source_ids,
    given_votable_get_df,
    given_source_ids_get_gaia_data,
    query_neighborhood,
    plot_group_neighborhood
)

from cdips.utils.get_vizier_catalogs import (
    get_k13_param_table,
    get_k13_member_stars
)

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

def main(overwrite=0):
    """
    first, try to acquire CG18 probabilities, by matching against the names
    of groups. (throws out PMS cands, and strings without any names).

    otherwise, acquire the K13 members (that were xmatched against Gaia DR2
    in the CDIPS I paper)

    otherwise, acquire the KC19 members

    acquire the mean properties of the group, and query the neighborhood
    based on those properties. the number of neighbor stars to randomly
    select is min(5* the number of group members, 5000). (cutoff group
    bounds based on parallax because further groups more uncertain).

    (note that this query produces hella field stars for the K13 members,
    beccause they're mostly contaminants!)
    """

    target_df = pd.DataFrame({
        'name':'',
        'toi_or_ticid':626.01,
        'source_id':5617241426979996800,
        'group_id':449
    }, index=[0])

    for ix, target_row in target_df.iterrows():

        ##########################################
        # acquire the CG18 info, if it is available.
        groupname = str(target_row['name']).split(',')[0]
        targetname = str(target_row['toi_or_ticid'])
        target_id = np.int64(target_row['source_id'])
        group_id = int(target_row['group_id'])

        group_in_cg18 = False
        group_in_k13 = False
        group_in_kc19 = True
        groupname = 'kc19group{:d}'.format(group_id)

        if not group_in_cg18 and not group_in_k13 and not group_in_kc19:
            print('found {} not in K13 or CG18 or KC19 groups. skip.'.
                  format(groupname))
            continue

        outpath = (
            '../results/neighborhoods/{}_{}_neighborhood.png'.
            format(targetname, groupname)
        )
        if not overwrite:
            if os.path.exists(outpath):
                print('found {}, skip'.format(outpath))
                continue

        if group_in_kc19:

            cutoff_probability = 1
            kc19_df = pd.read_csv('../data/kc19_string_table1.csv')
            kc19_df = kc19_df[kc19_df['group_id'] == group_id]
            group_source_ids = np.array(kc19_df['source_id']).astype(np.int64)

        else:
            raise NotImplementedError('this is a custom script for KELT25')


        ##########################################

        enforce_all_sourceids_viable = True

        group_df_dr2 = given_source_ids_get_gaia_data(
            group_source_ids, groupname, overwrite=overwrite,
            enforce_all_sourceids_viable=enforce_all_sourceids_viable
        )

        target_d = objectid_search(
            target_id,
            columns=('source_id', 'ra','dec', 'phot_g_mean_mag',
                     'phot_bp_mean_mag', 'phot_rp_mean_mag', 'l','b',
                     'parallax, parallax_error', 'pmra','pmra_error',
                     'pmdec','pmdec_error', 'radial_velocity'),
            forcefetch=True
        )
        target_df = pd.read_csv(target_d['result'])
        assert len(target_df) == 1

        # now acquire the mean properties of the group, and query the neighborhood
        # based on those properties. the number of neighbor stars to randomly
        # select is min(5* the number of group members, 5000). (cutoff group
        # bounds based on parallax because further groups more uncertain).
        bounds = {}
        params = ['parallax', 'ra', 'dec']

        plx_mean = group_df_dr2['parallax'].mean()
        if plx_mean > 5:
            n_std = 5
        elif plx_mean > 3:
            n_std = 4.5
        else:
            n_std = 4
        if group_in_kc19:
            n_std = 2
        if targetname == '1000.01':
            n_std = 5
        if targetname == '1097.01':
            n_std = 4
        print('bounding by {} stdevns'.format(n_std))

        for param in params:

            bounds[param+'_upper'] = (
                group_df_dr2[param].mean() + n_std*group_df_dr2[param].std()
            )

            bounds[param+'_lower'] = (
                group_df_dr2[param].mean() - n_std*group_df_dr2[param].std()
            )

        if bounds['parallax_lower'] < 0:
            bounds['parallax_lower'] = 0

        bounds['parallax_lower'] = 1.5

        assert bounds['ra_upper'] < 360
        assert bounds['ra_lower'] > 0
        assert bounds['parallax_lower'] >= 0

        n_max = min((50*len(group_df_dr2), 10000))
        nbhd_df = query_neighborhood(bounds, groupname, n_max=n_max,
                                     overwrite=overwrite,
                                     is_cg18_group=group_in_cg18,
                                     is_kc19_group=group_in_kc19,
                                     is_k13_group=group_in_k13)

        # ensure no overlap between the group members and the neighborhood sample.
        common = group_df_dr2.merge(nbhd_df, on='source_id', how='inner')
        snbhd_df = nbhd_df[~nbhd_df.source_id.isin(common.source_id)]

        n_std = 5
        if targetname == '1000.01':
            n_std = 8
        pmdec_min = group_df_dr2['pmdec'].mean() - n_std*group_df_dr2['pmdec'].std()
        pmdec_max = group_df_dr2['pmdec'].mean() + n_std*group_df_dr2['pmdec'].std()
        pmra_min = group_df_dr2['pmra'].mean() - n_std*group_df_dr2['pmra'].std()
        pmra_max = group_df_dr2['pmra'].mean() + n_std*group_df_dr2['pmra'].std()

        pmdec_min = min((pmdec_min, float(target_df['pmdec'])))
        pmdec_max = max((pmdec_max, float(target_df['pmdec'])))
        pmra_min = min((pmra_min, float(target_df['pmra'])))
        pmra_max = max((pmra_max, float(target_df['pmra'])))

        plot_group_neighborhood(targetname, groupname, group_df_dr2, target_df,
                                nbhd_df, cutoff_probability, extra_overplot=0,
                                pmdec_min=pmdec_min, pmdec_max=pmdec_max,
                                pmra_min=pmra_min, pmra_max=pmra_max,
                                group_in_k13=group_in_k13,
                                group_in_cg18=group_in_cg18,
                                group_in_kc19=group_in_kc19,
                                plx_ylim=[1.40,4.05])

        plot_group_neighborhood(targetname, groupname, group_df_dr2, target_df,
                                nbhd_df, cutoff_probability, extra_overplot=1,
                                pmdec_min=pmdec_min, pmdec_max=pmdec_max,
                                pmra_min=pmra_min, pmra_max=pmra_max,
                                group_in_k13=group_in_k13,
                                group_in_cg18=group_in_cg18,
                                group_in_kc19=group_in_kc19,
                                plx_ylim=[1.40,4.05])

if __name__ == "__main__":
    main()
