"""
One-time updates of stellar rotation parameters for planet candidates in
candidates.csv

The relevant manually updated parameters are:
    'rot_quality', 'Prot', 'vsini', 'rot_amp', 'Mp_pred'

update_single: Manually update 1 new candidate.
update_many: Manually update multiple new candidates.
"""

from cdips_followup.manage_candidates import update_candidate_rot_params
from cdips.utils import today_YYYYMMDD
import numpy as np, pandas as pd
import astropy.units as u

from astrobase.services.identifiers import gaiadr2_to_tic

def main():
    update_single = 1
    update_many = 0

    if update_single:
        # * rot_quality: 0 is good. 1 is so-so. 2 is very uncertain.
        # * ticid or source_id: str.
        # * At least one of Prot or vsini.
        #     (E.g., 3.14 and '--', or 42 and None, or 42 and 42)
        # * rot_amp: float
        # * Mp_pred (astropy quantity).

        update_candidate_rot_params(
            ticid='268301217',
            rot_quality='0',
            Prot=7.0*u.day,
            vsini=None,  # quantity, None, or '--'
            rot_amp=1.1e-2,
            Mp_pred=1*u.Mjup
        )

    if update_many:
        raise NotImplementedError

if __name__ == "__main__":
    main()
