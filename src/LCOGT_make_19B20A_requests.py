"""
Given targets, their ephemerides, and positions, submit requests to LCOGT to
get 1m and 2m imaging follow-up.
"""

###########
# imports #
###########
import requests
import numpy as np, pandas as pd

import datetime as dt
from astropy.time import Time
from astropy.coordinates import get_body, get_sun, get_moon, SkyCoord
import astropy.units as u

from astroplan import (FixedTarget, Observer, EclipsingSystem,
                       PrimaryEclipseConstraint, is_event_observable,
                       AtNightConstraint, AltitudeConstraint,
                       LocalTimeConstraint, MoonSeparationConstraint,
                       moon)

with open('/home/luke/.lcogt_api_token', 'r') as f:
    l = f.readlines()
token = str(l[0].replace('\n',''))

#############
# functions #
#############

def _given_Gmag_get_exptime_defocus(Gmag)

    #FIXME implement using the table

    return exptime, defocus

def make_request_group(targetname, ra, dec, Gmag, starttime, endtime,
                       max_airmass=2.5, min_lunar_distance=20, filtermode="ip",
                       telescope_class="1m0"):

    exptime, defocus = _given_Gmag_get_exptime_defocus(Gmag)

    API_TOKEN = token  # API token obtained from https://observe.lco.global/accounts/profile/
    PROPOSAL_ID = 'NOAO2019B-013'  # Proposal IDs may be found here: https://observe.lco.global/proposals/

    # starttime e.g., '2019-05-02 00:00:00',
    obsdate = (
          str(starttime).split('-')[0]      # year
        + str(starttime).split('-')[1]      # month
        + str(starttime).split('-')[2][:2]  # date
    )

    requestname = (
        '{targetname:s}_{telescope_class:s}_{obsdate:s}_{filtermode:s}_{exptime:d}_{defocus:.1f}'.
        format(targetname=targetname, telescope_class=telescope_class,
               obsdate=obsdate, filtermode=filtermode, exptime=exptime,
               defocus=defocus)
    )

    # The target of the observation
    target = {
        'name': str(targetname),
        'type': 'ICRS',
        'ra': float(ra), # decimal deg
        'dec': float(dec),
        'epoch': 2015.5
    }

    # Constraints used for scheduling the observation
    constraints = {
        'max_airmass': max_airmass,
        'min_lunar_distance': min_lunar_distance
    }

    # The configurations for this request. In this example we are taking 2 exposures with
    # different filters and exposure times. The fields acquisition_config and guiding_config 
    # are required fields in a configuration that are ultimately filled in with defaults 
    # if the submitted values are empty.
    configurations = [
        {
            'type': 'EXPOSE',
            'instrument_type': '1M0-SCICAM-SINISTRO',
            'target': target,
            'constraints': constraints,
            "instrument_configs": [
                {
                    "optical_elements": {
                        "filter": filtermode
                    },
                    "mode": "full_frame",
                    "exposure_time": int(exptime),
                    "exposure_count": int(expcount),
                    "bin_x": 1,
                    "bin_y": 1,
                    "rotator_mode": "",
                    "extra_params": {
                        "defocus": float(defocus)
                    }
                }
            ],
            'acquisition_config': {},
            'guiding_config': {
                "optional": False,
                "mode": "ON"
            },
            "priority": 1
        }
    ]

    # The time windows during which this request should be considered for observing. In this example
    # we only provide one. These times are in UTC.
    windows = [{
        'start': starttime # e.g., '2019-05-02 00:00:00',
        'end': endtime # '2019-05-30 00:00:00'
    }]

    # The telescope class that should be used for this observation
    location = {
        'telescope_class': telescope_class
    }

    # The full RequestGroup, with additional meta-data
    requestgroup = {
        'name': requestname,  # The title
        'proposal': PROPOSAL_ID,
        'ipp_value': 1.05,
        'operator': 'SINGLE',
        'observation_type': 'NORMAL',
        'requests': [{
            'configurations': configurations,
            'windows': windows,
            'location': location,
        }]
    }

    return requestgroup


def get_requests_given_ephem(targetname, ra, dec, Gmag,
                             period, period_unc, epoch, epoch_unc, depth,
                             depth_unc, duration, duration_unc,
                             min_search_time=Time(dt.datetime.today().isoformat()),
                             max_search_time=Time('2020-01-29 23:00:00'),
                             max_airmass=2.5, min_lunar_distance=20
                            ):
    """
    Given an ephemeris, and the basic details of a target, generate LCOGT
    requests for any available transits at the given site, between
    min_search_time and max_search_time.

    Allowed sites include Siding Spring Observator and CTIO.
    """

    # TODO: you need to make a loop here.
    make_request_group(targetname, ra, dec, Gmag, starttime, endtime,
                       max_airmass=2.5, min_lunar_distance=20, filtermode="ip",
                       telescope_class="1m0", max_airmass=max_airmass,
                       min_lunar_distance=min_lunar_distance)


def get_all_requests_19B():

    df = pd.read_csv('../data/20190912_19B20A_LCOGT_1m_2m.csv')
    #FIXME implement!!

    result = []

    for ix, r in df.iterrows():

        #FIXME IMPLEMENT
        this = get_requests_given_ephem(targetname, ra, dec, Gmag, period,
                                        period_unc, epoch, epoch_unc, depth,
                                        depth_unc, duration, duration_unc)

        result.append(this)

        import IPython; IPython.embed()

        assert 0


    return result


if __name__ == "__main__":

    r = get_all_requests_19B()
