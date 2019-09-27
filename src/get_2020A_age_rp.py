"""
use cdips env

for 2020A proposal figure, need to give best estimate at Rp and age for the
available systems.
"""
import os
from glob import glob
import pandas as pd, numpy as np
from cdips.utils import collect_cdips_lightcurves as ccl
from astropy import units as u

df = pd.read_csv('../results/20190912_2020A_targets_obs_info.csv', sep='|')

cdips_df = ccl.get_cdips_pub_catalog(ver=0.3)

#
# GET THE AGES
#
outpath = '../results/20190926_2020A_targets_age.csv'

if not os.path.exists(outpath):
    age = []
    for ix, r in df.iterrows():
        _age = r['reference_or_age']
        try:
            age.append(float(_age))
        except ValueError:
            mdf = cdips_df[cdips_df['source_id']==r['source_id']]

            if r['name']=='COL':
                age.append(np.log10(0.042e9)) # Gagne & Faherty 2018
                continue
            if r['name']=='PscEri':
                age.append(np.log10(0.120e9)) # Curtis+2019
                continue
            if 'Pozzo_1' in str(r['name']):
                age.append(np.log10(0.020e9)) # Franciosini+2018 == gamma vel
                continue
            if 'PMS' in str(r['reference_or_age']):
                age.append(np.log10(0.500e9)) # 500 Myr upper bound
                continue

            if len(mdf)==1:
                age.append(float(mdf['k13_logt']))
                continue
            else:
                age.append(np.nan)

    df['age'] = age
    df.to_csv(outpath, index=False)
else:
    df = pd.read_csv(outpath)

#
# NOW GET THE RADII
#

# uploaded converged objects
df1 = pd.read_csv(
    '/home/luke/Dropbox/proj/cdips/data/exoFOP_uploads/'
    '20190918_s6_and_s7_w_sourceid.csv', sep='|'
)

# TOIs
df2 = pd.read_csv(
    '/home/luke/Dropbox/proj/cdips/data/'
    'toi-plus-2019-08-29.csv'
)

# these from vetting reports
_d = {
5290781443841554432: (18.15, 3.00),
5618515825371166464: (14.99, 3.00),
5715454237977779968: (40.44, 8.00)
}
# all exofoptess, mostly eastman
_t = {
580.01: (11.77,1.79),
581.01: (4.71, 0.48),
861.01: ( (np.sqrt(4384e-6)*2.293*u.Rsun).to(u.Rearth).value,
          (0.4*np.sqrt(4384e-6)*2.293*u.Rsun).to(u.Rearth).value
        ),
588.01: (12.11,2.69),
1014.01: ( (np.sqrt(2800e-6)*2.1*u.Rsun).to(u.Rearth).value,
           (0.3*np.sqrt(2800e-6)*2.1*u.Rsun).to(u.Rearth).value
         ),
1034.01: ( (np.sqrt(510e-6)*1.644*u.Rsun).to(u.Rearth).value,
          (0.2*np.sqrt(510e-6)*1.644*u.Rsun).to(u.Rearth).value
         ),
1097.01: (2.35,0.12)
}

outpath = '../results/20190926_2020A_targets_age_rp.csv'

if not os.path.exists(outpath):
    rplanet = []
    rplanet_unc = []
    for ix, r in df.iterrows():

        mdf = df1[df1['source_id'] == r['source_id']]

        if len(mdf) == 0:

            m2df = df2[str(r['toi_or_ticid'])==df2['toi_id'].astype(str)]

            if len(m2df) == 1:
                if (np.isfinite(float(m2df['R_p'])) and
                    np.isfinite(float(m2df['R_p Err']))
                ):
                    rplanet.append(float(m2df['R_p']))
                    rplanet_unc.append(float(m2df['R_p Err']))
                    continue

            if str(r['source_id']) in np.array(list(_d.keys())).astype(str):
                rplanet.append(float(_d[r['source_id']][0]))
                rplanet_unc.append(float(_d[r['source_id']][1]))
                continue

            #print('{} in {}?'.format(str(r['toi_or_ticid']), np.array(list(_t.keys())).astype(str)))

            if str(r['toi_or_ticid']) in np.array(list(_t.keys())).astype(str):
                rplanet.append(float(_t[float(r['toi_or_ticid'])][0]))
                rplanet_unc.append(float(_t[float(r['toi_or_ticid'])][1]))
                continue

            rplanet.append(np.nan)
            rplanet_unc.append(np.nan)
            continue

        else:
            rplanet.append(float(mdf['radius']))
            rplanet_unc.append(float(mdf['radius_unc']))
            del mdf
            continue

    df['rplanet'] = rplanet
    df['rplanet_unc'] = rplanet_unc
    df.to_csv(outpath, index=False)

else:
    df = pd.read_csv(outpath)

