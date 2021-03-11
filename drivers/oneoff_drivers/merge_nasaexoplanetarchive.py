"""
Merge NASA exoplanet archive against CDIPS target list.
"""

from cdips.utils.tapqueries import get_exoplanetarchive_planetarysystems
from cdips.utils.catalogs import get_cdips_catalog
from cdips_followup.utils import get_cdips_candidates
from cdips_followup.paths import DATADIR, RESULTSDIR

import os
from glob import glob
import pandas as pd, numpy as np

#
# columns described at
# https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
#
df_ps = get_exoplanetarchive_planetarysystems(overwrite=0)

def get_sourceid(x):
    if len(x) > 0:
        return str(x.replace('Gaia DR2 ',''))
    else:
        return ''

source_id = df_ps.gaia_id.apply(get_sourceid)
df_ps['source_id'] = source_id

sdf_ps = df_ps.drop_duplicates('pl_name', keep='first')

df_cdips = get_cdips_catalog(ver=0.4)
df_cdips['source_id'] = df_cdips.source_id.astype(str)

mdf = sdf_ps.merge(df_cdips, how='left', on='source_id')

assert len(mdf) == len(sdf_ps)

outdf = mdf[~pd.isnull(mdf.reference)]

outpath = os.path.join(DATADIR, '20210311_ps_X_cdipsv0pt4.csv')
outdf.to_csv(outpath, index=False)

selcols = ['pl_name', 'discoverymethod', 'source_id', 'phot_g_mean_mag', 'cluster', 'reference',
           'pl_rade', 'pl_orbper', 'pl_masse']
sdf = outdf[selcols].sort_values(by='pl_orbper')
sdf = sdf[sdf.reference != 'Zari_2018_UMS']
outpath = os.path.join(RESULTSDIR,
                       '20210311_NASA_ExoplanetArchive_ps_X_cdips0pt4',
                       '20210311_ps_X_cdipsv0pt4_short.csv')
sdf.to_csv(outpath, index=False)
