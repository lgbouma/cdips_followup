from astropy.time import Time
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.LCOGT_dedicated_requests import (
    get_dedicated_request,
    given_dedicated_requests_validate_submit
)

if __name__ == "__main__":

    ####################
    overwrite = 1
    validate = 1
    submit = 0
    max_duration_error = 20
    savstr = '20191210_2m_requests' # eg, 20191207_TOI1098_request_2m_tc_secondary. "ephemupdate" if it is one. (this cancels pending observations)
    source_id = '2919143383943171200'
    manual_ephemeris = False
    max_search_time = Time('2020-03-01 23:59:00')
    filtermode = 'ip'
    telescope_class = '2m0'

    # submit_eventclasses = ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']
    # submit_eventclasses = ['OIBEO']
    submit_eventclasses = ['OIBEO']
    # submit_eventclasses = ['OIB', 'BEO']

    create_eventclasses = ['OIBEO'] #, 'IBEO', 'OIBE', 'OIB', 'BEO']
    # create_eventclasses = ['OIBEO', 'IBEO', 'OIBE']
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

    given_dedicated_requests_validate_submit(
        requests, submit_eventclasses, validate=validate, submit=submit,
        max_duration_error=max_duration_error)
