"""
Merge PATHOS-II against CDIPS candidates.
"""
from cdips_followup.utils import get_cdips_candidates
import pandas as pd
import numpy as np
from astrobase.services.identifiers import tic_to_gaiadr2
import os
from glob import glob

pathos_df = pd.read_csv('../../data/Nardiello_2020_PATHOS-II_candidates.csv')
cdips_df = get_cdips_candidates()

mdf = pathos_df.merge(cdips_df, how='inner', left_on='TICID', right_on='ticid')

print(f'N_pathos: {len(pathos_df)}')
print(f'N_inner: {len(mdf)}')

# PATHOS object ids that I don't already have a classification on
# (see 20200610_pathos_merge_notes.txt)
objs = np.array([3, 8, 9, 15, 23, 31])

n_int_path = '../../data/Nardiello_2020_PATHOS-II_interesting.csv'
if not os.path.exists(n_int_path):

    sdf = pathos_df[pathos_df.PATHOSID.isin(objs)]

    gaia_ids = []

    for t in sdf.TICID:

        g = tic_to_gaiadr2(t)
        gaia_ids.append(g)

    sdf['GAIA'] = gaia_ids
    sdf.to_csv(n_int_path, index=False)
else:
    sdf = pd.read_csv(n_int_path)


classfxdir = '/Users/luke/Dropbox/proj/cdips/results/vetting_classifications'
classpaths = glob(os.path.join(classfxdir, '*LGB*classifications.csv'))

cdf = pd.concat([pd.read_csv(f) for f in classpaths])

gaia_ids = [n.split('_')[4].split('-')[0].replace('gaiatwo','').lstrip('0') for
           n in cdf.Name]
cdf['GAIA'] = np.array(gaia_ids).astype(np.int64)

smdf = sdf.merge(cdf, how='left', on='GAIA')

import IPython; IPython.embed()
