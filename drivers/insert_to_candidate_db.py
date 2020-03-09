"""
One-time inserts of new young planet candidate to candidates.csv

insert_single: Manually insert 1 new candidate.
insert_many: Manually insert multiple new candidates.
"""
from cdips_followup.manage_candidates import insert_candidate
from cdips.utils import today_YYYYMMDD
import numpy as np, pandas as pd

from astrobase.services.identifiers import gaiadr2_to_tic

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
        'init_priority': init_priority,
        'current_priority': current_priority,
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


def manual_insert_many():

    sourceids = [
        '5514577414951631488',
        '1835201042675810688',
        '2192733400508791552'
    ]
    comments = [
        "PC. T=9, Rp=4.6, P=3.5, Teff 8500K.",
        "KC19, needs check. (Interesting cluster). FIES (Buchave) showed vrot",
        "200km/s. T=8, Rp=4.4, P=2.3d.",
        "High priority, VPC.",
        "TRES FR, weak odd-even, Rp=5, P=1.1",
        "Hot star, needs phot FU"
    ]

    for s, c in zip(sourceids, comments):
        ticid = gaiadr2_to_tic(s)

        insert_single_candidate(
            ticid = ticid,
            comment = c,
            nbhd_rating = 0,
            init_priority = 0,
            pending_photometry_observations = 'to_assess',
            pending_spectroscopic_observations = 'to_assess'
        )


if __name__ == "__main__":

    insert_single = 0
    insert_many = 1

    if insert_single:

        insert_single_candidate(
            ticid = '62483237',
            comment = 'Whopping Ca HK.',
            nbhd_rating = 0,
            init_priority = 0,
            pending_photometry_observations = 'SG1 has lost',
            pending_spectroscopic_observations = 'Veloce17'
        )

    if insert_many:
        manual_insert_many()
