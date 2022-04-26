"""
go from the 20190912_toi_youngstar_xmatching.2020A_targets spreadsheet to
something with sufficient information to do some followup.
"""
import os
import numpy as np, matplotlib.pyplot as plt, pandas as pd

from astropy.io.votable import from_table, writeto, parse
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

from astroquery.mast import Catalogs
from astroquery.gaia import Gaia
from astrobase.services.gaia import objectid_search

# updated at 20190919
df = pd.read_csv('../data/20190912_toi_youngstar_xmatching - 2020A_targets.csv')

columns=('source_id', 'ra','dec', 'parallax', 'phot_g_mean_mag',
         'phot_bp_mean_mag', 'phot_rp_mean_mag')

gaia_dfs, tic_dfs = [], []

for ix, r in df.iterrows():

    source_id = r['source_id']

    # get the gaia information
    gaia_d = objectid_search(source_id, columns=columns)
    gaia_df = pd.read_csv(gaia_d['result'])

    ra, dec = float(gaia_df['ra']), float(gaia_df['dec'])
    targetcoord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    radius = 1.0*u.arcminute

    stars = Catalogs.query_region(
        "{} {}".format(float(targetcoord.ra.value), float(targetcoord.dec.value)),
        catalog="TIC",
        radius=radius
    )

    # require finite dr2 value
    sel = ~stars['GAIA'].mask
    stars = stars[sel]

    # TICv8: just get whichever star matches the source_id!
    target_star = stars[stars['GAIA']==str(source_id)]

    if len(target_star) != 1:
        raise AssertionError('didnt get match for {}'.format(source_id))

    cols =  ['ID', 'Bmag', 'Vmag', 'Jmag', 'Hmag', 'Kmag', 'Tmag', 'Teff', 'logg', 'rad', 'mass']

    tic_df = target_star[cols].to_pandas()

    gaia_dfs.append(gaia_df)
    tic_dfs.append(tic_df)

gdf = pd.concat(gaia_dfs)
tdf = pd.concat(tic_dfs)
tdf = tdf.rename(columns={"ID": "ticid", "Teff": "ticteff", "rad":"ticrstar"})
tdf['source_id'] = gdf['source_id']

outdf = df.merge(gdf, on='source_id', how='left')
outdf = outdf.merge(tdf, on='source_id', how='left')

outpath = '../results/20190912_2020A_targets_obs_info.csv'
outdf.to_csv(outpath, index=False, sep='|')
print('made {}'.format(outpath))
