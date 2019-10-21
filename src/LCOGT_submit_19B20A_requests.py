"""
submit requests made in LCOGT_make_19B20A_requests.py
"""
##########
import pickle, requests, socket

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
##########

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

    billed_durn = (
        requestgroup_dict['request_durations']['requests'][0]['duration']
    )

    start = Time(requestgroup['requests'][0]['windows'][0]['start'])
    end = Time(requestgroup['requests'][0]['windows'][0]['end'])
    window_durn = (end - start).value*24*60*60

    if (window_durn - billed_durn)/60 > max_duration_error:

        errmsg = (
            'ERROR! got a window of {:.2f} min; billed {:.2f} min.'.
            format(window_durn/60, billed_durn/60)
        )

        raise AssertionError(errmsg)

    else:

        return 1



def submit_single_request(requestgroup):

    if validate_single_request(requestgroup):
        pass

    # Submit the fully formed RequestGroup
    response = requests.post(
        'https://observe.lco.global/api/requestgroups/',
        headers={'Authorization': 'Token {}'.format(token)},
        json=requestgroup  # Make sure you use json!
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


def submit_all_requests(max_N_transit_per_object=2):

    pkl_savpath = (
        '../results/LCOGT_19B20A_observability/all_requests_19B.pkl'
    )
    with open(pkl_savpath, 'rb') as f:
        r = pickle.load(f)

    df = pd.read_csv(
        '../results/LCOGT_19B20A_observability/all_requests_summary.csv'
    )

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

    for _r, ind in zip(r, time_sort_inds):

        _requests_sorted = np.array(_r)[ind]

        _requests_to_submit = _requests_sorted[:max_N_transit_per_object]

        for requestgroup in _requests_to_submit:

            print(requestgroup)

            # FIXME change to submit
            # submit_single_request(requestgroup)




if __name__=="__main__":
    submit_all_requests()
