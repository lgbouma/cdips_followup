"""
Given a source_id, make LCOGT photometry followup requests, and optionally
submit them to the LCOGT API.
"""
import numpy as np
from astropy.time import Time
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.LCOGT_dedicated_requests import (
    get_dedicated_request,
    given_dedicated_requests_validate_submit
)

def main():

    ####################
    savstr = '20200208_1m_requests' # eg, 20191207_TOI1098_request_2m_tc_secondary. "ephemupdate" if it is one. (this cancels pending observations)
    overwrite = 1
    validate = 1
    submit = 1

    source_id = '5519619186857962112' # TOI580
    max_n_events = 10 # else None. n_events is per eventclass.

    filtermode = 'ip'
    telescope_class = '1m0'

    create_eventclasses = ['OIBEO']
    # create_eventclasses = ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']
    # create_eventclasses = ['OIBEO', 'IBEO', 'OIBE']
    # create_eventclasses = ['OIB', 'BEO']

    submit_eventclasses = ['OIBEO']
    # submit_eventclasses = ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']
    # submit_eventclasses = ['OIBEO', 'IBEO', 'OIBE']
    # submit_eventclasses = ['OIB', 'BEO']

    max_duration_error = 20
    max_search_time = Time('2020-05-29 23:59:00')
    manual_ephemeris = False

    ####################

    if manual_ephemeris:
        period = 42
        epoch = 2458660.00000
        duration = 2.00000

    else:
        # get ephemeris from ephemerides.csv
        d = query_ephemeris(source_id=source_id)
        period, epoch, duration = d['period'], d['epoch'], d['duration']

    # "requests" is a list of lists. Higher level is each eventclass. Level
    # below is each event, in that eventclass.
    requests = get_dedicated_request(savstr, source_id, period, epoch,
                                     duration, create_eventclasses,
                                     overwrite=overwrite,
                                     max_search_time=max_search_time,
                                     filtermode=filtermode,
                                     telescope_class=telescope_class)

    # if a maximum number of events is set, impose it!
    if isinstance(max_n_events, int):

        _requests = []
        for ix in range(len(create_eventclasses)):
            print('starting with {} {} events.'.
                  format(len(requests[ix]), create_eventclasses[ix])
            )

        for eventclass in requests:
            _eventclass = []
            starttimes = []
            for req in eventclass:
                starttimes.append(req['requests'][0]['windows'][0]['start'])

            # sort by start time, cut to get the closest ones.
            sort_times = np.sort(starttimes)
            sel_times = sort_times[ : max_n_events]

            for req in eventclass:
                starttime = req['requests'][0]['windows'][0]['start']
                if starttime in sel_times:
                    _eventclass.append(req)

            if len(_eventclass) > 0:
                _requests.append(_eventclass)

        assert len(_requests[0]) <= max_n_events
        requests = _requests

        print('WRN!: trimmed to {} events.'.format(len(requests[0])))
        if len(sel_times)>0:
            print('WRN!: max time: \n{}'.format(repr(sel_times[-1])))
            print('\nWRN!: selected times: \n{}'.format(repr(sel_times)))
        else:
            print('WRN!: got no times')


    given_dedicated_requests_validate_submit(
        requests, submit_eventclasses, validate=validate, submit=submit,
        max_duration_error=max_duration_error)


if __name__ == "__main__":
    main()
