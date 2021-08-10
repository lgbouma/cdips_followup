"""
Merge NASA exoplanet archive, or the TOI list, against the CDIPS target list.
"""

from cdips.utils.catalogs import (
    get_cdips_catalog, get_exofop_toi_catalog, TODAYSTR
)
from cdips_followup.paths import DATADIR, RESULTSDIR
from astrobase.services.identifiers import tic_to_gaiadr2

import os
from glob import glob
import pandas as pd, numpy as np

OVERWRITE = 0
CDIPSVER = 0.6
verstr = '2021-08-07'

#
# latest ExoFOP TOI catalog
#
savstr = verstr.replace('-','')
outdir = os.path.join(RESULTSDIR, 'KNOWNPLANET_YOUNGSTAR_XMATCH',
                      f'{savstr}_TOI_x_CDIPStargets')
if not os.path.exists(outdir):
    os.mkdir(outdir)
outpath  = os.path.join(outdir, 'exofop_CTOI_cache.csv')
if not os.path.exists(outpath):
    df = get_exofop_toi_catalog(ver=verstr)
    #... this runs ~15 TICID/minute. 4 hours. annoying.
    # alternative ideas:
    # 1. run a smart MAST query, uploading a table of TIC IDs, xmatching on the
    # Gaia ID column of the TIC8
    # 2. bulk download all TIC8. 
    source_ids = []
    for ix, r in df.iterrows():
        print(ix, r['TOI'])
        try:
            source_ids.append(tic_to_gaiadr2(str(r['TIC ID'])))
        except:
            source_ids.append(np.nan)
    df['source_id'] = source_ids
    df.to_csv(outpath, index=False)
df = pd.read_csv(outpath)

df['source_id'] = df.source_id.astype(str)

df_cdips = get_cdips_catalog(ver=CDIPSVER)
df_cdips['source_id'] = df_cdips.source_id.astype(str)

mdf = df.merge(df_cdips, how='left', on='source_id')

outdf = mdf[~pd.isnull(mdf.reference_id)]

outpath = os.path.join(
    outdir, f'{savstr}_ps_X_cdipsv{str(CDIPSVER).replace(".","pt")}.csv'
)
outdf.to_csv(outpath, index=False)

selcols = ['source_id', 'TIC ID', 'TOI', 'phot_g_mean_mag',
           'cluster', 'reference_id', 'mean_age', 'Planet Radius (R_Earth)',
           'Period (days)', 'Depth (ppm)', 'Comments']
sdf = outdf[selcols].sort_values(by='TOI')

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
    &
    (sdf.reference_id != 'Zari2018ums')
    &
    (sdf.reference_id != 'Oh2017')
)

_sel = (sdf.reference_id.str.contains('HATSandHATNcandidates20210505'))
_df = sdf[_sel]
outpath = os.path.join(
    outdir, f'{savstr}_ps_X_cdipsv{str(CDIPSVER).replace(".","pt")}_HATSHATN.csv')
_df.to_csv(outpath, index=False)
print(f'Made {outpath}')

sdf = sdf[sel]
outpath = os.path.join(
    outdir, f'{savstr}_ps_X_cdipsv{str(CDIPSVER).replace(".","pt")}_short.csv')
sdf.to_csv(outpath, index=False)
print(f'Made {outpath}')

sel &= (sdf.reference_id != 'HATSandHATNcandidates20210505')
sdf = sdf[sel]
outpath = os.path.join(
    outdir, f'{savstr}_ps_X_cdipsv{str(CDIPSVER).replace(".","pt")}_shorter.csv')
sdf.to_csv(outpath, index=False)
print(f'Made {outpath}')
