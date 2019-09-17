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

def main(overwrite=1):

    target_df = pd.read_csv('../data/20190912_youngstar_cands_with_gaia.csv')

    # first, try to acquire CG18 probabilities, by matching against the names
    # of groups. (throws out PMS cands, and strings without any names).
    starget_df = target_df[~pd.isnull(target_df['name'])]

    Vizier.ROW_LIMIT = 3000
    cg18_vizier_str ='J/A+A/618/A93'
    catalog_list = Vizier.find_catalogs(cg18_vizier_str)
    catalogs = Vizier.get_catalogs(catalog_list.keys())

    cg18_tab1 = catalogs[0]
    cg18_tab1_df = cg18_tab1.to_pandas()

    cg18_tab1_df['Cluster'] = cg18_tab1_df['Cluster'].str.decode('utf-8')
    cg18_clusters = np.array(cg18_tab1_df['Cluster']).astype(str)

    for ix, target_row in target_df.iterrows():

        ##########################################
        # acquire the CG18 info, if it is available. (if not, skip this row,
        # for now..)
        groupname = str(target_row['name']).split(',')[0]
        if groupname not in cg18_clusters:
            print('found {} not in CG18 clusters. skip.'.format(groupname))
            continue
        outpath = '../results/{}_neighborhood.png'.format(groupname)
        #if os.path.exists(outpath):
        #    print('found {}, skip'.format(outpath))
        #    continue

        targetname = str(target_row['toi_or_ticid'])
        target_id = np.int64(target_row['source_id'])

        cutoff_probability = 0.1

        # now get the CG18 members for this group.
        v = Vizier(column_filters={"Cluster":groupname})
        v.ROW_LIMIT = 3000
        catalog_list = v.find_catalogs(cg18_vizier_str)
        catalogs = v.get_catalogs(catalog_list.keys())

        group_tab = catalogs[1]
        group_df = group_tab.to_pandas()
        assert len(group_df) < v.ROW_LIMIT
        np.testing.assert_array_equal(group_tab['Source'], group_df['Source'])

        group_source_ids = np.array(group_df['Source']).astype(np.int64)

        group_df = group_df[group_df['PMemb'] > cutoff_probability]
        ##########################################

        group_df_dr2 = given_source_ids_get_gaia_data(group_source_ids, groupname,
                                                      overwrite=overwrite)

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

        n_max = min((50*len(group_df_dr2), 2000))
        nbhd_df = query_neighborhood(bounds, groupname, n_max=n_max,
                                     overwrite=overwrite, is_cg18_group=True,
                                     is_kc19_group=False)

        # ensure no overlap between the group members and the neighborhood sample.
        common = group_df_dr2.merge(nbhd_df, on='source_id', how='inner')
        snbhd_df = nbhd_df[~nbhd_df.source_id.isin(common.source_id)]

        pmdec_min = group_df_dr2['pmdec'].mean() - 2*group_df_dr2['pmdec'].std()
        pmdec_max = group_df_dr2['pmdec'].mean() + 2*group_df_dr2['pmdec'].std()
        pmra_min = group_df_dr2['pmra'].mean() - 2*group_df_dr2['pmra'].std()
        pmra_max = group_df_dr2['pmra'].mean() + 2*group_df_dr2['pmra'].std()

        pmdec_min = min((pmdec_min, target_df['pmdec']))
        pmdec_max = max((pmdec_max, target_df['pmdec']))
        pmra_min = min((pmdec_min, target_df['pmra']))
        pmra_max = max((pmdec_max, target_df['pmra']))

        plot_group_neighborhood(targetname, groupname, group_df_dr2, target_df,
                                nbhd_df, cutoff_probability, extra_overplot=0,
                                pmdec_min=pmdec_min, pmdec_max=pmdec_max,
                                pmra_min=pmra_min, pmra_max=pmra_max)

        plot_group_neighborhood(targetname, groupname, group_df_dr2, target_df,
                                nbhd_df, cutoff_probability, extra_overplot=1,
                                pmdec_min=pmdec_min, pmdec_max=pmdec_max,
                                pmra_min=pmra_min, pmra_max=pmra_max)

if __name__ == "__main__":
    main()
