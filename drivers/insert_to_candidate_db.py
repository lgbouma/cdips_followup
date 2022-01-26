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
        #'2890062263458259968',
        #'2072431332214824576',
        #'1833519030401513984',
        '5404435280772350848'
    ]

    comments = [
        #"V-shaped. Cool if planet (probably is not, since V-shaped).  Neighborhood plot looks fine -- Platais5 is definitely a real thing though Prot cleanup required.   Star Prot 2.5 or 3 days from SPOC LC.  Teff 3250K, it's an M dwarf.  Smearing expected here on transit, given the 2 hour transit?  Needs some RVs to verify.  SG1 finds VPC, detecting event on target.",
        #"Good youngest HJ cand.  Teff 7800K; an early F/late A-star.  Rp=22Re.  One could set a mass upper limit via RVs, and get a DT observation to confirm.  Cluster membership from Gaia blind kinematics is good (actual cluster plot would be better).  Needs recon spectroscopy to check if high RV scatter."
        #"This candidate HJ orbits a Teff 5850K G-dwarf.  The star's rotation period appears to be 2.5-3 days. i.e., consistent with an age of ~100 Myr give or take.  And the transit looks really good.  Kinematically, the str is in the outskirts of the Alessi12 cluster, at least as discussed by KC19.  It is one of the closest stars reported in the cluster, and so its proper motions are pretty different from the other cluster members.  Needs specFU (e.g., for Li) to verify youth."
        "Probably multi. Prot 3d? is consistent with the youth, b/c Teff 6500K... but also consistent with a good deal older than 100 Myr.   Neighborhood plot looks ok.  Definitely worth getting Li, spectra, etc to check"
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

    insert_single = 1
    insert_many = 0

    if insert_single:

        insert_single_candidate(
            ticid = '56551765',
            comment = 'V1096 Tau = Herbig Anon 1.  Taurus 1-3 Myr.  Maybe HEB system, maybe 1-3Myr HJ..  Need to check for chromaticity, then careful RVs to untangle the SB2. ',
            nbhd_rating = 0,
            init_priority = 0,
            pending_photometry_observations = 'PP0.',
            pending_spectroscopic_observations = 'SP1.'
        )

    if insert_many:
        manual_insert_many()
