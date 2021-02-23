"""
LCOGT_dedicated_requests.py
Luke Bouma, 2018-2020

Contents:
    get_dedicated_request
    given_dedicated_requests_validate_submit
"""
import pandas as pd, numpy as np
import datetime as dt
from astropy.time import Time
import pickle, os
from copy import deepcopy

from cdips_followup.LCOGT_make_requests import (
    make_single_request_from_row,
    get_current_semester_str
)

from cdips_followup.LCOGT_submit_requests import (
    validate_single_request,
    submit_single_request
)

from cdips_followup.manage_ephemerides import query_ephemeris

from cdips_followup import __path__
DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data')
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')

def get_dedicated_request(savstr, source_id, period, epoch, duration,
                          eventclasses, overwrite=0, semesterstr=None,
                          min_search_time=Time(dt.datetime.today().isoformat()),
                          max_search_time=None, filtermode='ip',
                          telescope_class='1m0',
                          sites=None
                         ):
    """
    #
    # savstr: e.g., request_2m_tc_secondary. "ephemupdate" if it is one...
    #
    """

    assert isinstance(source_id, str)

    if semesterstr is None:
        semesterstr = get_current_semester_str()

    r = {'source_id': source_id, 'period': period,
         'epoch': epoch, 'duration': duration}

    init_savstr = f'{savstr}_all_classes'
    _savstr = deepcopy(savstr)

    if not 'ephemupdate' in init_savstr:
        resultsdir = (
            os.path.join(RESULTSDIR, f'LCOGT_{semesterstr}_observability/')
        )
    else:
        resultsdir = (
            os.path.join(RESULTSDIR, f'LCOGT_{semesterstr}_updated_requests/')
        )
    if not os.path.exists(resultsdir):
        os.mkdir(resultsdir)

    pkl_savpath = (
        os.path.join(resultsdir, f'{init_savstr}.pkl')
    )

    if not overwrite and os.path.exists(pkl_savpath):
        with open(pkl_savpath, 'rb') as f:
            requests = pickle.load(f)
    else:
        requests = []
        for eventclass in eventclasses:
            savstr = '{}_{}'.format(_savstr, eventclass)
            df = pd.DataFrame(r, index=[0])
            for _, row in df.iterrows():
                req = make_single_request_from_row(
                    row, savstr, eventclass, min_search_time=min_search_time,
                    max_search_time=max_search_time, filtermode=filtermode,
                    telescope_class=telescope_class, sites=sites
                )
            requests.append(req)

        with open(pkl_savpath, 'wb') as f:
            pickle.dump(requests, f, pickle.HIGHEST_PROTOCOL)
            print('saved {:s}'.format(pkl_savpath))

    return requests


def given_dedicated_requests_validate_submit(requests,
                                             submit_eventclasses=None,
                                             validate=1, submit=0,
                                             max_duration_error=15,
                                             overwrite_acceptability=None,
                                             overwrite_ipp=None,
                                             raise_error=True):
    """
    Choose desired events to submit, then validate and submit them.

    submit_eventclasses : list
        example: ['OIBEO', 'BEO', 'IBEO']

    Assumes name pattern following "TIC308538095.01_OIBEO_1m0_20200220_ip_90_2.0"
    """

    requestgroups = []
    for r in requests:
        if len(r) > 0:
            for _r in r:
                this_name = _r['name']
                this_eventclass = this_name.split('_')[1]
                if this_eventclass in submit_eventclasses:
                    requestgroups.append(_r)

    for requestgroup in requestgroups:

        if isinstance(overwrite_acceptability, int):
            requestgroup['acceptability_threshold'] = overwrite_acceptability
        if isinstance(overwrite_ipp, int):
            requestgroup['ipp_value'] = overwrite_ipp

        if validate:
            if not submit:
                print(requestgroup)
            requestgroup, is_modified = (
                validate_single_request(
                    requestgroup, max_duration_error=max_duration_error,
                    raise_error=raise_error
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

        if submit:
            if isinstance(requestgroup, dict):
                print('SUBMITTING...')
                print(requestgroup)
                submit_single_request(requestgroup)
            else:
                print('vvv DID NOT SUBMIT B/C FAILED TO VALIDATE vvv')
                print(requestgroup)
                print('^^^ DID NOT SUBMIT B/C FAILED TO VALIDATE ^^^')
