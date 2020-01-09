"""
One-time initialization of candidates.csv
"""
from cdips_followup.manage_candidates import insert_candidate
import os, shutil
from glob import glob
import numpy as np, pandas as pd

homedir = os.path.expanduser('~')

src = '../data/candidate_database/BASETEMPLATE.csv'
dst = '../data/candidate_database/candidates.csv'
shutil.copyfile(src, dst)
print('WRN! REBASED {}'.format(dst))

#
# import the s6+s7, + ~fall 2019 TOI xmatch candidates as exported from the
# deprecated google spreadsheet
#
df = pd.read_csv(
    '../data/candidate_database/CDIPS_followup_tracker_20200108_import_s6-s7_plus_TOI.csv'
)
df.source_id = df.source_id.astype(str)

for ix, r in df.iterrows():

    d = {
        'nbhd_rating': r.nbhd_rating,
        'init_priority': r.init_priority,
        'current_priority': r.current_priority,
        'pending_spectroscopic_observations': r.pending_spectroscopic_observations,
        'pending_photometry_observations': r.pending_photometry_observations,
        'comment': r.comment,
        'candidate_provenance': '20191110_followup_comment_strings',
        'isretired': r.isretired
    }

    insert_candidate(source_id=r.source_id, manual_dict=d)

#
# import the s8-s11 candidates. for priorities, use the average from the
# classification files
#
df = pd.read_csv(
    '../data/candidate_database/'
    'CDIPS_followup_tracker_20200109_import_s8-s11.csv'
)
df.source_id = df.source_id.astype(str)

classglob = os.path.join(
    homedir,
    'Dropbox/proj/cdips/results',
    'vetting_classifications/sector-*_PCs_CLEAR_THRESHOLD.csv'
)
priority_df = pd.concat(
    [pd.read_csv(f, sep=';') for f in glob(classglob)]
)
priority_df['source_id'] = list(map(
    lambda x: str(x.split('gaiatwo')[-1].split('-')[0].lstrip('0')),
    list(priority_df.Name)
))

priority_df_no_dups = priority_df.groupby('source_id').mean().reset_index()

mdf = df.merge(priority_df_no_dups, on='source_id', how='left')
assert len(mdf) == len(df)

mdf['init_priority'] = 2 - mdf.average_score
mdf['current_priority'] = 2 - mdf.average_score

for ix, r in mdf.iterrows():

    d = {
        'nbhd_rating': r.nbhd_rating,
        'init_priority': r.init_priority,
        'current_priority': r.current_priority,
        'pending_spectroscopic_observations': r.pending_spectroscopic_observations,
        'pending_photometry_observations': r.pending_photometry_observations,
        'comment': r.comment,
        'candidate_provenance': 'CDIPS_followup_tracker_20200109_import_s8-s11',
        'isretired': 0
    }

    #FIXME FIXME: YOU NEED TO DO SOMETHING WHEN YOU HAVE MATCHING SOURCE_IDS.
    # LIKE RAISE AN ERROR (BY DEFAULT)

    insert_candidate(source_id=r.source_id, manual_dict=d)
