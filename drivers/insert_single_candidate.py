"""
One-time inserts of new young planet candidate to candidates.csv
"""
from cdips_followup.manage_candidates import insert_candidate
from cdips.utils import today_YYYYMMDD
import numpy as np, pandas as pd

def insert_single_candidate(
    ticid = None,
    comment = '',
    nbhd_rating = 1,
    init_priority = 1,
    pending_photometry_observations = '--',
    pending_spectroscopic_observations = '--'):


    current_priority = init_priority

    d = {
        'nbhd_rating': 1,
        'init_priority': 1,
        'current_priority': 1,
        'pending_spectroscopic_observations': pending_spectroscopic_observations,
        'pending_photometry_observations': pending_photometry_observations,
        'comment': comment,
        'candidate_provenance': today_YYYYMMDD()+'_manual_insert',
        'isretired': 0
    }

    insert_candidate(
        source_id=None, ticid=ticid, manual_dict=d,
        raise_error_if_duplicate=True
    )


if __name__ == "__main__":

    insert_single_candidate(
        ticid = 'ENTER TICID',
        comment = 'MAKE A COMMENT',
        nbhd_rating = 0,
        init_priority = 1,
        pending_photometry_observations = '--',
        pending_spectroscopic_observations = '--'
    )
