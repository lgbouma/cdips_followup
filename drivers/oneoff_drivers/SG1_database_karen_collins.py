"""
For SG1 followup, Karen Collins requires 

TIC
NumberPlanet # (start at 0.11, 0.12, ...)
Pipeline
Discovery Sector
Priority (<4 Re = 2) (>4 Re = 3)"Re
Disposition (starts at "PC")PC"
T0 (BTJD) (Reference Epoch)"Epoch"
T0 Uncertainty (Days)"Days
Period (Days)Days"Period Uncertainty"Uncertainty
Duration (hrs)hrs
Depth (ppm)ppm
Rp (R_Earth)R_Earth"
Vetting Comments and/or Observational Needs"Needs"
RA (Decimal Degrees)"Degrees
Dec (Decimial Degrees)Degrees
Depth Uncertainty (ppm)ppm
Duration Uncertainty (hrs)
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
#
# 467522800.01 through 1030389830.01, taken from the candidate db
#
source_ids = [5240531632175135616, 5246508676134913024, 5252320149182351744,
              5254794943702272640, 5255504399239636992, 5256502966294604800,
              5256717749007641344, 5304507540851982592, 5308445029440522112,
              5312011948237411456, 5312141072137904512, 5321581959996214016,
              5322083917107299712, 5329853517263005696, 5338645555866749952,
              5339389268061191040, 5404579488593432576, 5405530566146401408,
              5406188928778889472, 5489726768531119616, 5514373043228306816,
              5523449717870971776, 5544257185155788672, 5548219137873425280,
              5790098845432414464, 5838183443852841216, 5869166341556812032 ]

##########################################
##########################################

#
# s8-s11 candidates from chat_database_rafael.py
#
dfpath = (
    '../../results/SG1_candidate_sharing/20200121_s8_to_s11_cdips_cands.csv'
)
df = pd.read_csv(dfpath, sep='|')
df['source_id'] = df.source_id.astype(str)

scols = [
    'source_id', 'ticid', 'epoch', 'epoch_unc', 'period', 'period_unc',
    'duration', 'duration_unc', 'rp', 'rp_unc', 'comment', 'depth', 'depth_unc'
]
kdf = df[scols]

ras, decs = [], []
for ix, source_id in enumerate(np.array(kdf.source_id)):
    print('{}/{}'.format(ix, len(kdf)))
    ra,dec = given_sourceid_get_radec(source_id)
    ras.append(ra)
    decs.append(dec)

kdf['ra'] = ras
kdf['dec'] = decs

kdf['NumberPlanet'] = '0.11'
kdf['Pipeline'] = 'CDIPS'
kdf['DiscoverySector'] = 'S8-S11' #FIXME: do this via tesscut last sector
kdf['Priority'] = 3
kdf['Disposition'] = 'PC'

colorder = [
    'ticid', 'NumberPlanet', 'Pipeline', 'DiscoverySector', 'Priority',
    'Disposition', 'epoch', 'epoch_unc', 'period', 'period_unc', 'duration',
    'duration_unc', 'rp', 'rp_unc', 'comment', 'ra', 'dec', 'depth',
    'depth_unc'
]

kdf = kdf[colorder]

# epoch in BJTD = BJD - 2457000
kdf['epoch'] -= 2457000

outpath = (
    '../../results/SG1_candidate_sharing/20200201_s8_to_s11_cdips_cands_karen_collins.csv'
)
kdf.to_csv(outpath, index=False, sep='|')
print('made {}'.format(outpath))
