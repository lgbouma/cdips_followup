"""
Validate and submit requests made in LCOGT_make_19B20A_requests.py
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

if socket.gethostname() == 'brik':
    api_file = '/home/luke/.lcogt_api_token'
elif 'astro' in socket.gethostname():
    api_file = '/Users/luke/.lcogt_api_token'
else:
    raise NotImplementedError('where to get API file?')

with open(api_file, 'r') as f:
    l = f.readlines()
token = str(l[0].replace('\n',''))

#############
# functions #
#############

def validate_single_request(requestgroup, max_duration_error=15):
    """
    Submit the RequestGroup through the "validate" API, cf.
    https://developers.lco.global/#validate-a-requestgroup

    max_duration_error: in minutes, is the maximum allowable difference between
    the start & end times of the request, and the _billed duration_ of the
    request. By design in the API, the billed duration is always shorter than
    the (end-start) time. I allotted 1 hour on either side for scheduling, so a
    bit of slack on either is fine.
    """

    is_modified = False

    response = requests.post(
        'https://observe.lco.global/api/requestgroups/validate/',
        headers={'Authorization': 'Token {}'.format(token)},
        json=requestgroup
    )

    # Make sure the API call was successful
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        print('API call failed: {}'.format(response.content))
        raise exc

    requestgroup_dict = response.json()

    # If you get an error because your incorrectly estimated the number of
    # exposures, correct it here.
    if len(requestgroup_dict['errors']) >= 1:

        if 'non_field_errors' in requestgroup_dict['errors']:

            print(42*'-')
            print('GOT ERROR: {}'.
                  format(requestgroup_dict['errors']['non_field_errors']))
            print(42*'-')

            return np.nan, np.nan

        errmsg = (
            requestgroup_dict['errors']['requests'][0]['non_field_errors'][0]
        )

        if 'the target is visible for a maximum of' in errmsg:

            # get the strings of durations, and decrement the requested number
            # of exposures by the right multiple!
            sr = search("According{}maximum of {} hours "
                        "within{}your request {} hours. Consider{}",
                        errmsg)

            max_dur = float(sr[1])
            req_dur = float(sr[3])

            if req_dur == max_dur:
                # {:.1f} formatted strings. genius ._.
                req_dur += 0.01

            if not req_dur > max_dur:
                errmsg = (
                    'ERR! max dur: {}, req dur: {}'.format(max_dur, req_dur)
                )
                raise ValueError(errmsg)

            diff_dur_sec = (req_dur - max_dur)*60*60

            # previously, guessed
            #
            # expcount = np.floor(
            #     (endtime-starttime).to(u.hr)
            #     /
            #     (exptime*u.second + read_time_per_exposure).to(u.hr)
            # )
            #
            # that produced the difference above...
            exptime_sec = (
                requestgroup['requests'][0]['configurations'][0]['instrument_configs'][0]['exposure_time']
            )

            expcount = (
                requestgroup['requests'][0]['configurations'][0]['instrument_configs'][0]['exposure_count']
            )

            read_time_per_exposure = 30*u.second # from Bayliss' completed runs
            n_exposures_diff = int(
                np.ceil(diff_dur_sec/
                        (exptime_sec + read_time_per_exposure.value)
                )
            )

            new_expcount = expcount - n_exposures_diff

            print(42*'-')
            print('WRN!: max durn: {} hr, req durn: {} hr. had {} exposures, decrement to {}'.
                  format(max_dur, req_dur, expcount, new_expcount))
            print(42*'-')
            requestgroup['requests'][0]['configurations'][0]['instrument_configs'][0]['exposure_count'] = new_expcount

            is_modified = True

            return requestgroup, is_modified

        else:
            raise NotImplementedError('got new API error: {}'.format(errmsg))

    billed_durn = (
        requestgroup_dict['request_durations']['requests'][0]['duration']
    )

    start = Time(requestgroup['requests'][0]['windows'][0]['start'])
    end = Time(requestgroup['requests'][0]['windows'][0]['end'])
    window_durn = (end - start).value*24*60*60

    expcount = (
        requestgroup['requests'][0]['configurations'][0]['instrument_configs'][0]['exposure_count']
    )

    if (window_durn - billed_durn)/60 > max_duration_error:

        errmsg = (
            'ERROR! got a window of {:.2f} min; billed {:.2f} min.'.
            format(window_durn/60, billed_durn/60)
        )

        print(42*'-')
        print(errmsg)
        print(42*'-')
        #import IPython; IPython.embed()
        #raise AssertionError(errmsg) #FIXME
        return np.nan, np.nan

    else:

        print(42*'-')
        print('ACCEPTED! window durn: {:.2f} min, billed {:.2f} min. had {:d} exposures'.
              format(window_durn/60, billed_durn/60, expcount))
        print(42*'-')
        return requestgroup, is_modified


def submit_single_request(requestgroup):

    # Submit the fully formed RequestGroup
    response = requests.post(
        'https://observe.lco.global/api/requestgroups/',
        headers={'Authorization': 'Token {}'.format(token)},
        json=requestgroup
    )

    # Make sure the API call was successful
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        print('API call failed: {}'.format(response.content))
        raise exc

    # The API returns the newly submitted requestgroup as json
    requestgroup_dict = response.json()

    # Print out the url on the portal where we can view the submitted request
    print('View the observing request: '
          'https://observe.lco.global/requestgroups/{}/'.
          format(requestgroup_dict['id']))


def submit_all_requests(savstr, validate_all=1, submit_all=0,
                        max_N_transit_per_object=3, max_duration_error=15):
    """
    savstr: used for directory management

    validate_all: if true, first validates observation requests to ensure that
    they can be submitted

    submit_all: actually submits them

    max_N_transit_per_object:

    max_duration_error: in minutes, maximum acceptable difference between
    _desired_ observation window, and the window that the LCOGT system accepts.
    """

    if submit_all:
        assert validate_all

    if not 'ephemupdate' in savstr:
        resultsdir = '../results/LCOGT_19B20A_observability/'
    else:
        resultsdir = '../results/LCOGT_19B20A_updated_requests/'

    pkl_savpath = (
        os.path.join(resultsdir, '{}.pkl'.format(savstr))
    )
    mult_savpath = (
        os.path.join(resultsdir, '{}_summary.csv'.format(savstr))
    )

    with open(pkl_savpath, 'rb') as f:
        r = pickle.load(f)

    df = pd.read_csv(mult_savpath)
    if submit_all:
        print('ATTEMPTING TO SUBMIT THE FOLLOWING')

        df['submit_durn'] = (
            df['sched_duration'] *
            np.minimum(df['n_requests'], max_N_transit_per_object)
        )

        print(df)

        print(42*'=')
        print('\nTotal time: {:.1f} hr\n'.format(np.sum(df['submit_durn'])))
        print(42*'=')

    #
    # sort all the available transit windows for each target by time. submit
    # the earliest `max_N_transit_per_object' (e.g., 2 transits).
    #
    starts = []
    for _r in r:
        starts.append(
            [ Time(__r['requests'][0]['windows'][0]['start']) for __r in _r ]
        )

    time_sort_inds = []
    for start in starts:
        time_sort_inds.append(
            np.argsort(start)
        )

    #
    # iterate over available requests for each target that met the magnitude
    # and depth cuts 
    #
    for _r, ind in zip(r, time_sort_inds):

        _requests_sorted = np.array(_r)[ind]

        _requests_to_submit = _requests_sorted[:max_N_transit_per_object]

        for requestgroup in _requests_to_submit:

            if '451.01' in requestgroup['name']: #FIXME FIXME
                print('WRN! SKIPPING 451.01 BECAUSE ITS A WEIRDO')
                continue
            if 'TIC53682439.01' in requestgroup['name']: #FIXME FIXME
                print('WRN! SKIPPING 5368* BECAUSE ITS A WEIRDO')
                continue


            if validate_all:
                if not submit_all:
                    print(requestgroup)
                requestgroup, is_modified = (
                    validate_single_request(
                        requestgroup, max_duration_error=max_duration_error
                    )
                )

                n_iter = 0
                if is_modified and np.isfinite(is_modified):
                    while is_modified:
                        if n_iter >= 10:
                            raise AssertionError('too many iterations')
                        requestgroup, is_modified = (
                            validate_single_request(
                                requestgroup,
                                max_duration_error=max_duration_error
                            )
                        )

                        if not isinstance(requestgroup, dict):
                            if not np.isfinite(requestgroup):
                                break

                        n_iter += 1

            if submit_all:
                if isinstance(requestgroup, dict):
                    print('SUBMITTING...')
                    print(requestgroup)
                    submit_single_request(requestgroup)
                else:
                    print('vvv DID NOT SUBMIT B/C FAILED TO VALIDATE vvv')
                    print(requestgroup)
                    print('^^^ DID NOT SUBMIT B/C FAILED TO VALIDATE ^^^')


if __name__=="__main__":

    validate_all = 1
    submit_all = 1
    max_N_transit_per_object = 2
    max_duration_error = 20

    eventclass = 'OIBEO'
    savstr = 'request_19B_59859387_{}'.format(eventclass)

    # eventclass = 'OIB'
    # savstr = 'request_19B_2m_faint_{}'.format(eventclass)

    # eventclass = 'OIBE'
    # savstr = 'bright_shallow_19B_{}'.format(eventclass)

    # eventclass = 'OIB'
    # savstr = 'midpartials_19B_{}'.format(eventclass)

    # eventclass = 'OIB'
    # savstr = 'toppartials_19B_{}'.format(eventclass)

    # max_duration_error = 15
    # savstr = 'request_TIC29786532_19B'
    # max_N_transit_per_object = 2

    # savstr = 'request_19B_2m_faint_v2'
    # max_N_transit_per_object = 2

    # savstr = 'request_19B_2m_faint'
    # max_N_transit_per_object = 4 # actually 3, b/c one fails

    # savstr = 'all_requests_19B_easyones'
    # max_N_transit_per_object = 3

    submit_all_requests(savstr, validate_all=validate_all,
                        submit_all=submit_all,
                        max_N_transit_per_object=max_N_transit_per_object,
                        max_duration_error=max_duration_error)
