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

OVERWRITE = 0
CDIPSVER = 0.5

#
# columns described at
# https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
#
df_ps = get_exoplanetarchive_planetarysystems(overwrite=OVERWRITE)

def get_sourceid(x):
    if len(x) > 0:
        return str(x.replace('Gaia DR2 ',''))
    else:
        return ''

source_id = df_ps.gaia_id.apply(get_sourceid)
df_ps['source_id'] = source_id

sdf_ps = df_ps.drop_duplicates('pl_name', keep='first')

df_cdips = get_cdips_catalog(ver=CDIPSVER)
df_cdips['source_id'] = df_cdips.source_id.astype(str)

mdf = sdf_ps.merge(df_cdips, how='left', on='source_id')

assert len(mdf) == len(sdf_ps)

outdf = mdf[~pd.isnull(mdf.reference_id)]

outpath = os.path.join(DATADIR, '20210510_ps_X_cdipsv0pt5.csv')
outdf.to_csv(outpath, index=False)

selcols = ['pl_name', 'discoverymethod', 'source_id', 'phot_g_mean_mag',
           'cluster', 'reference_id', 'mean_age', 'pl_rade', 'pl_orbper', 'pl_masse']
sdf = outdf[selcols].sort_values(by='pl_orbper')

ref_ids = np.array(sdf.reference_id)
from collections import Counter
res = Counter(ref_ids)
print(res.most_common(n=20))

sel = (
    (sdf.reference_id != 'NASAExoArchive_ps_20210506')
    &
    (sdf.reference_id != 'HATSandHATNcandidates20210505,NASAExoArchive_ps_20210506')
    &
    (sdf.reference_id != 'NASAExoArchive_ps_20210506,HATSandHATNcandidates20210505')
    &
    (sdf.reference_id != 'Zari2018ums,NASAExoArchive_ps_20210506')
    &
    (sdf.reference_id != 'NASAExoArchive_ps_20210506,Zari2018ums')
    &
    (sdf.reference_id != 'Oh2017,NASAExoArchive_ps_20210506')
    &
    (sdf.reference_id != 'NASAExoArchive_ps_20210506,Oh2017')
)

sdf = sdf[sel]
outdir = os.path.join(RESULTSDIR,
                       '20210510_NASA_ExoplanetArchive_ps_X_cdips0pt5')
if not os.path.exists(outdir):
    os.mkdir(outdir)
outpath = os.path.join(outdir, '20210510_ps_X_cdipsv0pt5_short.csv')
sdf.to_csv(outpath, index=False)
print(f'Made {outpath}')
