"""
To insert CDIPS candidates into the "CHIRON planning" spreadsheet,
a particular format is required...

TICID
Rp P Teff TIC RA DEC Vmag N_exp
"""

##########
# config #
##########

import numpy as np, pandas as pd
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.utils import (
    get_cdips_candidates,
    given_sourceid_get_radec
)
from astroquery.gaia import Gaia
from numpy import array as nparr
from astropy.coordinates import ICRS
from astropy import units as u

df = get_cdips_candidates()

# select * where K<2 and not L contains 'SP2' and not L contains '--'
sel = (
    (df.current_priority < 2)
    &
    ~(df.pending_spectroscopic_observations.str.contains('SP2'))
    &
    ~(df.pending_spectroscopic_observations.str.contains('--'))
)
df = df[sel]

scols = ['source_id', 'ticid', 'targetid', 'rp', 'period', 'tic_Vmag',
         'tic_teff', 'pending_spectroscopic_observations']
df = df[scols]

ras, decs = [], []
for ix, source_id in enumerate(nparr(df.source_id.astype(str))):
    print('{}/{}'.format(ix, len(df)))
    ra, dec = given_sourceid_get_radec(source_id)
    ras.append(ra)
    decs.append(dec)

c = ICRS(nparr(ras)*u.deg, nparr(decs)*u.deg)

rahmsstr = c.ra.to_string(u.hour, sep=' ', pad=True)
decdmsstr = c.dec.to_string(u.degree, sep=' ', pad=True)

df['ra'] = rahmsstr
df['dec'] = decdmsstr

# Rp P Teff TIC RA DEC Vmag N_exp
df['N_exp'] = 3 # to start with

df['Notes'] = 'CDIPS. ' + df.pending_spectroscopic_observations
df = df.drop(['pending_spectroscopic_observations'], axis=1)

outpath = '../../results/SG1_candidate_sharing/20200122_quinn_chiron_export.csv'
df.to_csv(outpath, index=False, sep='|')
print('made {}'.format(outpath))
