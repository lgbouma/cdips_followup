import pandas as pd, numpy as np
import os
from cdips_followup.manage_candidates import save_candidates_csv_file

from cdips.utils.catalogs import (
    get_toi_catalog,
    get_exofop_ctoi_catalog
)


df = pd.read_csv('../../data/candidate_database/candidates.csv', sep="|")
df.targetid = df.targetid.astype(str)

#
# get rp uncertainties that I reported to exofop-tess
#
ctoidf = get_exofop_ctoi_catalog()
ctoidf = ctoidf[(ctoidf.User == 'bouma')]
ctoidf.CTOI = ctoidf.CTOI.astype(str)

mdf_ctoi = df.merge(ctoidf, how='left', left_on='targetid', right_on='CTOI')
assert len(mdf_ctoi) == len(df)

rp_ctoi = np.array(mdf_ctoi['Radius (R_Earth)'])
rp_unc_ctoi = np.array(mdf_ctoi['Radius (R_Earth) Error'])

#
# get rp uncertainties from the MIT-QLP TEV list
#
toidf = get_toi_catalog()
toidf['Full TOI ID'] = toidf['Full TOI ID'].astype(str)
mdf_toi = df.merge(toidf, how='left', left_on='targetid',
                   right_on='Full TOI ID')
assert len(mdf_toi) == len(df)

rp_toi = np.array(mdf_toi['Planet Radius Value'])
rp_unc_toi = np.array(mdf_toi['Planet Radius Error'])

#
# add the radii uncertainties to the candidates database. 
# TEV takes priority.
# ensure that rp is nan iff. rp_unc is also nan.
#
rp_unc_to_db = -1*np.ones_like(rp_ctoi)
rp_unc_to_db[~pd.isnull(rp_unc_toi)] = rp_unc_toi[~pd.isnull(rp_unc_toi)]
rp_unc_to_db[~pd.isnull(rp_unc_ctoi)] = rp_unc_ctoi[~pd.isnull(rp_unc_ctoi)]

df['rp_unc'] = rp_unc_to_db

np.testing.assert_array_equal(df.rp == -1, df.rp_unc==-1)

save_candidates_csv_file(df)
