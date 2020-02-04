"""
Given updated ephemerides for a target (e.g., after getting new data!), cancel
whatever requests were pending, recalculate the request windows, and resubmit
them.
"""

###########
# imports #
###########
import pickle, requests, socket
from parse import search

import os
from glob import glob

import numpy as np, pandas as pd

from astropy.time import Time
import astropy.units as u

from cdips.utils import today_YYYYMMDD

from cdips_followup.LCOGT_make_requests import make_all_request_files
from cdips_followup.LCOGT_submit_requests import submit_all_requests

if socket.gethostname() == 'brik':
    api_file = '/home/luke/.lcogt_api_token'
elif 'astro' in socket.gethostname():
    api_file = '/Users/luke/.lcogt_api_token'
else:
    raise NotImplementedError('where to get API file?')

with open(api_file, 'r') as f:
    l = f.readlines()
token = str(l[0].replace('\n',''))

from cdips_followup import __path__
DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data')
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')

#############
# functions #
#############

def get_requests(limit=500, user='lbouma'):
    # limit: maximum number of requests returned for pagination purposes
    # identifier: TICID, or TOI id. (Leading string)

    response = requests.get(
        'https://observe.lco.global/api/requestgroups/?limit={}&state=PENDING&user={}'.
        format(limit, user),
        headers={'Authorization': 'Token {}'.format(token)}
    )

    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        print('API call failed: {}'.format(response.content))
        raise exc

    response_dict = response.json()

    if not response_dict['count'] == len(response_dict['results']):
        wrn = (
            'WRN! found {} requests, but pagination only collected {}'.
            format(response_dict['count'],
                   len(response_dict['results']))
        )
        print(wrn)

    return response_dict


def cancel_single_request(requestid):

    response = requests.post(
        'https://observe.lco.global/api/requestgroups/{}/cancel/'.
        format(requestid),
        headers={'Authorization': 'Token {}'.format(token)}
    )

    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        print('API call failed: {}'.format(response.content))
        raise exc

    response_dict = response.json()

    return response_dict


def get_new_ephemeris(targetid, verbose=True):
    # read in the new ephemeris provided by Joel Hartman
    ephem_path = (
        os.path.join(DATADIR,
                     'updated_ephemerides/{}/{}.updateephem.txt'.
                     format(today_YYYYMMDD(), targetid))
    )

    with open(ephem_path, 'r') as f:
        lines = f.readlines()

    epoch = [l for l in lines if '- Epoch: ' in l][0]
    period = [l for l in lines if '- Period: ' in l][0]
    dur = [l for l in lines if '- Transit duration: ' in l][0]

    if verbose:
        print(epoch, period, dur)

    epoch = float(search('{} - Epoch: {} +/- {}', epoch)[1].strip())
    period = float(search('{} - Period: {} +/- {}', period)[1].strip())
    dur = float(search('{} - Transit duration: {} +/- {}', dur)[1].strip())

    if verbose:
        print(epoch, period, dur)

    ephem_dict = {
        'period':period,
        'epoch':epoch,
        'duration':dur
    }

    return ephem_dict


def get_requestids_given_targetid(requestdict, targetid):

    requestids = []

    for r in requestdict['results']:

        if str(targetid) in r['name']:
            requestids.append(r['id'])

    return requestids


def cancel_pending_requests(targetid):
    # cancel pending requests for a single target, e.g., "TIC53682439.01"

    rd = get_requests()

    to_cancel_ids = get_requestids_given_targetid(rd, targetid)

    #
    # cancel all pending requests for the targets to update
    #
    if len(to_cancel_ids)>0:
        print('CANCELLING ALL PENDING REQUESTS FOR {}'.format(targetid))
        for to_cancel_id in to_cancel_ids:
            cancel_single_request(to_cancel_id)

    else:
        print('Did not find pending requests for {}'.format(targetid))


def prepare_for_ephem_update(targets_to_update, target_is_faint=False):
    # cancel pending requests. make all potential request dictionaries.

    make_eventclasses = ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']

    rd = get_requests()

    for targetid in targets_to_update:

        ephem_dict = get_new_ephemeris(targetid)

        to_cancel_ids = get_requestids_given_targetid(rd, targetid)

        #
        # cancel all pending requests for the targets to update
        #
        if len(to_cancel_ids)>0:
            print('CANCELLING ALL PENDING REQUESTS FOR {}'.format(targetid))
            for to_cancel_id in to_cancel_ids:
                cancel_single_request(to_cancel_id)

        else:
            print('Did not find pending requests for {}'.format(targetid))


        #
        # make all potential request dictionaries to submit
        #
        for eventclass in make_eventclasses:

            # note: the targetid as the first entry is used in "get_targets" to
            # find auxiliary information about the target.
            savstr = (
                '{}_ephemupdate_{}_on{}_19B'.
                format(targetid, eventclass, today_YYYYMMDD())
            )
            if target_is_faint:
                savstr = (
                    '{}_ephemupdate_2m_{}_on{}_19B'.
                    format(targetid, eventclass, today_YYYYMMDD())
                )

            outdir = (
                os.path.join(RESULTSDIR,
                             'LCOGT_19B20A_updated_requests/{}'.
                             format(savstr))
            )

            if not os.path.exists(outdir):
                _ = make_all_request_files(
                    savstr=savstr, overwrite=0, eventclass=eventclass,
                    ephem_dict=ephem_dict
                )
            else:
                print('found {}, continue'.format(outdir))



def submit_ephem_update(targetid, validate_all=None, submit_all=None,
                        eventclass=None, max_N_transit_per_object=2,
                        target_is_faint=False):

    max_duration_error = 20

    savstr = (
        '{}_ephemupdate_{}_on{}_19B'.
        format(targetid, eventclass, today_YYYYMMDD())
    )

    if target_is_faint:
        savstr = (
            '{}_ephemupdate_2m_{}_on{}_19B'.
            format(targetid, eventclass, today_YYYYMMDD())
        )

    submit_all_requests(
        savstr,
        validate_all=validate_all,
        submit_all=submit_all,
        max_N_transit_per_object=max_N_transit_per_object,
        max_duration_error=max_duration_error,
        target_is_faint=target_is_faint
    )


def main():

    prepare = 0
    validate = 1
    submit = 1
    target_is_faint = 0 # set to 1 if you want to do LCOGT/2m requests

    submit_eventclasses = ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']
    max_N_transit_per_object = 2

    # 20191116
    targets_to_prepare = [
        "TIC53682439"
        "TIC59859387"
    ]
    targets_to_submit = [
        "TIC53682439"
        "TIC59859387"
    ]

    ##########################################
    # # 20191108
    # targets_to_prepare = [
    #     "TIC238611475"
    # ]

    # # 20191030
    # targets_to_prepare = [
    #     "TIC308538095",
    #     "TIC59859387",
    #     "TIC349118653"
    # ]
    ##########################################

    if prepare:
        prepare_for_ephem_update(targets_to_prepare,
                                 target_is_faint=target_is_faint)

    if validate or submit:

        for targetid in targets_to_submit:

            for eventclass in submit_eventclasses:

                submit_ephem_update(
                    targetid,
                    validate_all=validate,
                    submit_all=submit,
                    eventclass=eventclass,
                    max_N_transit_per_object=max_N_transit_per_object,
                    target_is_faint=target_is_faint
                )


if __name__ == "__main__":

    main()
