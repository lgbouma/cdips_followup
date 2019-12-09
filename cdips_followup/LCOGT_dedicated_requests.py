import pandas as pd, numpy as np
import pickle, os
from copy import deepcopy

from cdips_followup.LCOGT_make_requests import make_single_request_from_row

from cdips_followup.LCOGT_submit_requests import (
    validate_single_request,
    submit_single_request
)

from cdips_followup.manage_ephemerides import query_ephemeris

def get_dedicated_request(savstr, source_id, period, epoch, duration,
                          eventclasses, overwrite=0, semesterstr='20A',
                          max_search_time=None):
    #
    # savstr: e.g., request_2m_tc_secondary. "ephemupdate" if it is one...
    #

    assert isinstance(source_id, str)

    r = {'source_id': source_id, 'period': period,
         'epoch': epoch, 'duration': duration}

    init_savstr = '{}_all_classes'.format(savstr)
    _savstr = deepcopy(savstr)

    if not 'ephemupdate' in init_savstr:
        resultsdir = (
            '../results/LCOGT_{}_observability/'.format(semesterstr)
        )
    else:
        resultsdir = (
            '../results/LCOGT_{}_updated_requests/'.format(semesterstr)
        )

    pkl_savpath = (
        os.path.join(resultsdir, '{}.pkl'.format(init_savstr))
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
                    row, savstr, eventclass, max_search_time=max_search_time
                )
            requests.append(req)

        with open(pkl_savpath, 'wb') as f:
            pickle.dump(requests, f, pickle.HIGHEST_PROTOCOL)
            print('saved {:s}'.format(pkl_savpath))

    return requests

def given_dedicated_requests_validate_submit(requests, validate=1, submit=0,
                                             max_duration_error=15,
                                             overwrite_acceptability=None,
                                             overwrite_ipp=None):

    requestgroups = []
    for r in requests:
        if len(r) > 0:
            for _r in r:
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

        if submit:
            if isinstance(requestgroup, dict):
                print('SUBMITTING...')
                print(requestgroup)
                submit_single_request(requestgroup)
            else:
                print('vvv DID NOT SUBMIT B/C FAILED TO VALIDATE vvv')
                print(requestgroup)
                print('^^^ DID NOT SUBMIT B/C FAILED TO VALIDATE ^^^')
