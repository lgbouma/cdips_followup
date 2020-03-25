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
    update_single = 0
    update_many = 1

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
            Mp_pred=1*u.Mjup,
            Tdur=1.38*u.hour
        )

    if update_many:

        #
        # this list was a one-off manual update from SP0 targets on 2020/03/24
        #
        ticids = ['268301217', '257605131', '460205581', '360630575',
                  '59859387', '383390264', '166527623', '62483237',
                  '217933560', '198262850']
        rot_qualitys = ['0', '0', '0', '0',
                        '1', '0', '0', '0',
                        '0', '0']
        Prots = [7*u.day, 4.5*u.day, None, 6.5*u.day,
                 None, 1.8*u.day, 1.5*u.day, 8*u.day,
                 3*u.day, 2.2*u.day]
        vsinis = [8*u.km/u.s, None, 17*u.km/u.s, None,
                  27.2*u.km/u.s, None, None, None,
                  None, None]
        rot_amps = [1.50E-02, 2.00E-02, 2.00E-02, 1.50E-02,
                    1.00E-03, 5.00E-03, 2.00E-02, 4.00E-03,
                    8.00E-02, 4.10E-02]
        Mp_preds = [1*u.Mjup, 15.6*u.Mearth, 49.6*u.Mearth, 6*u.Mearth,
                    20*u.Mearth, 9.5*u.Mearth, 20*u.Mearth, 10*u.Mearth,
                    1*u.Mjup, 15*u.Mearth]
        Tdurs = [1.38*u.hour, 4.436*u.hour, 1.775*u.hour, 2.343*u.hour,
                 1.73*u.hour, 2.64*u.hour, 3.1*u.hour, 2.136*u.hour,
                 2.54*u.hour, 1.92*u.hour ]

        for ticid, rot_quality, Prot, vsini, rot_amp, Mp_pred, Tdur in zip(
            ticids, rot_qualitys, Prots, vsinis, rot_amps, Mp_preds, Tdurs
        ):
            update_candidate_rot_params(
                ticid=ticid,
                rot_quality=rot_quality,
                Prot=Prot,
                vsini=vsini,
                rot_amp=rot_amp,
                Mp_pred=Mp_pred,
                Tdur=Tdur
            )

if __name__ == "__main__":
    main()
