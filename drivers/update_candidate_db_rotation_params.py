"""
One-time updates of stellar rotation parameters for planet candidates in
candidates.csv. Used for RM predictions.

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
            ticid='359357695',
            rot_quality='1',
            Prot=3.0*u.day,#None
            vsini=None,#9*u.km/u.s,  # quantity, None, or '--'
            rot_amp=0.02,
            Mp_pred=16*u.Mearth,#*u.Mearth,
            Tdur=1.2*u.hour
        )

    if update_many:

        #
        # this list was a one-off manual update from SP0 targets on 2020/03/24
        #
        ticids = [
            '180987952', '130415266', '154671430',
            '462162963', '140978827', '238597883',
            '117789567', '67600790'
        ]
        rot_qualitys = [
            '0', '1', '0',
            '0', '0', '1',
            '1', '0'
        ]
        Prots = [
            0.5*u.day, 0.5*u.day, 5*u.day,
            1*u.day, None, None,
            None, 1*u.day
        ]
        vsinis = [
            200*u.km/u.s, 200*u.km/u.s, 12*u.km/u.s,
            50*u.km/u.s, 13*u.km/u.s, 100*u.km/u.s,
            200*u.km/u.s, None
        ]
        rot_amps = [
            0.001, 0.001, 0.001,
            0.001, 0.002, 0.001,
            0.001, 0.05
        ]
        Mp_preds = [
            20*u.Mearth, 317.83*u.Mearth, 17*u.Mearth,
            317.83*u.Mearth, 317.83*u.Mearth, 20*u.Mearth,
            317.83*u.Mearth, 30*u.Mearth
        ]
        Tdurs = [
            1.514*u.hour, 6.1*u.hour, 2*u.hour,
            1.67*u.hour, 2.03*u.hour, 2.57*u.hour,
            1.33*u.hour, 1.22*u.hour
        ]

        # ticids = ['268301217', '257605131', '460205581', '360630575',
        #           '59859387', '383390264', '166527623', '62483237',
        #           '217933560', '198262850']
        # rot_qualitys = ['0', '0', '0', '0',
        #                 '1', '0', '0', '0',
        #                 '0', '0']
        # Prots = [7*u.day, 4.5*u.day, None, 6.5*u.day,
        #          None, 1.8*u.day, 1.5*u.day, 8*u.day,
        #          3*u.day, 2.2*u.day]
        # vsinis = [8*u.km/u.s, None, 17*u.km/u.s, None,
        #           27.2*u.km/u.s, None, None, None,
        #           None, None]
        # rot_amps = [1.50E-02, 2.00E-02, 2.00E-02, 1.50E-02,
        #             1.00E-03, 5.00E-03, 2.00E-02, 4.00E-03,
        #             8.00E-02, 4.10E-02]
        # Mp_preds = [1*u.Mjup, 15.6*u.Mearth, 49.6*u.Mearth, 6*u.Mearth,
        #             20*u.Mearth, 9.5*u.Mearth, 20*u.Mearth, 10*u.Mearth,
        #             1*u.Mjup, 15*u.Mearth]
        # Tdurs = [1.38*u.hour, 4.436*u.hour, 1.775*u.hour, 2.343*u.hour,
        #          1.73*u.hour, 2.64*u.hour, 3.1*u.hour, 2.136*u.hour,
        #          2.54*u.hour, 1.92*u.hour ]

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
