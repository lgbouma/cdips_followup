"""
For targets identified in the 2020A CHIRON and PFS proposals, we should observe
first with LCOGT. Most of the TOIs are too hot for PFS, but plausible for
CHIRON.
"""

import numpy as np, pandas as pd

tois = np.array([
    451.01,
    580.01,
    581.01,
    837.01,
    861.01,
    520.01,
    588.01,
    1014.01,
    1034.01,
    1097.01
]).astype(float)

df = pd.read_csv('../data/toi-plus-2019-10-19.csv')

_df = pd.DataFrame({'toi_id':tois})
mdf = df.merge(_df, how='inner', on='toi_id')

cols = [
    'tic_id', 'toi_id', 'Period', 'Period Err', 'Epoc', 'Epoc Err', 'Transit Depth',
    'Transit Depth Err', 'Duration', 'Duration Err'
]

odf = mdf[cols].sort_values(by='toi_id')

outpath = '../results/19B20A_LCOGT_1m_2m_ephemerides_TOIs-ONLY.csv'
odf.to_csv(outpath, index=False)
print('made {}'.format(outpath))
