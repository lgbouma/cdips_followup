"""
Merge stars with Kepler data against the CDIPS target list.
"""

from cdips.utils.catalogs import get_cdips_catalog
from cdips_followup.paths import DATADIR, RESULTSDIR

import os
from glob import glob
import pandas as pd, numpy as np
from astropy.table import Table

OVERWRITE = 0
CDIPSVER = 0.6

# Bedell's 1 arcsecond DR2 X Kepler match
KEPGAIAFUN = '/Users/luke/local/kepler_gws/kepler_dr2_1arcsec.fits'
df_kep = Table.read(KEPGAIAFUN, format='fits').to_pandas()
df_kep.source_id = df_kep.source_id.astype(str)
df_kep['planet?'] = df_kep['planet?'].str.decode('utf-8')

df_cdips = get_cdips_catalog(ver=CDIPSVER)
df_cdips['source_id'] = df_cdips.source_id.astype(str)

mdf = df_kep.merge(df_cdips, how='inner', on='source_id')

outdf = mdf[~pd.isnull(mdf.reference_id)]

outpath = os.path.join(DATADIR, '20210916_kepgaiafun_X_cdipsv0pt6.csv')
outdf.to_csv(outpath, index=False)
print(f'Wrote {outpath}')

selcols = ['kepid', 'source_id', 'kepmag', 'phot_g_mean_mag_x',
           'cluster', 'reference_id', 'mean_age', 'planet?', 'nkoi', 'nconfp',
           'kepler_gaia_ang_dist']

sdf = outdf[selcols].sort_values(by='mean_age')
sdf = sdf.rename({'phot_g_mean_mag_x':'phot_g_mean_mag'},axis='columns')

ref_ids = np.array(sdf.reference_id)
from collections import Counter
res = Counter(ref_ids)
print(42*'-')
print('all kepler stars...')
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
                       '20210916_kepgaiafun_X_cdips0pt6')
if not os.path.exists(outdir):
    os.mkdir(outdir)
outpath = os.path.join(outdir, '20210916_kepgaiafun_X_cdipsv0pt6_short.csv')
sdf.to_csv(outpath, index=False)
print(f'Made {outpath}')


sdf0 = sdf[(sdf.nkoi > 0)]

ref_ids = np.array(sdf0.reference_id)
res = Counter(ref_ids)
print(42*'-')
print('NKOI > 0')
print(res.most_common(n=20))

outpath = os.path.join(outdir, '20210916_kepgaiafun_X_cdipsv0pt6_KOIs_only.csv')
sdf0.to_csv(outpath, index=False)
print(f'Made {outpath}')

