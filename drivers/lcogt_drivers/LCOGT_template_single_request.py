"""
Given a source_id, make LCOGT photometry followup requests, and optionally
submit them to the LCOGT API.
"""
from astropy.time import Time
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.LCOGT_dedicated_requests import (
    get_dedicated_request,
    given_dedicated_requests_validate_submit
)

def main():

    ####################
    overwrite = 1
    validate = 1
    submit = 0

    savstr = '20200115_1m_requests' # eg, 20191207_TOI1098_request_2m_tc_secondary. "ephemupdate" if it is one. (this cancels pending observations)
    source_id = '5240531632175135616'

    filtermode = 'ip'
    telescope_class = '1m0'
    max_duration_error = 20
    max_search_time = Time('2020-04-15 23:59:00')
    manual_ephemeris = False

    create_eventclasses = ['OIBEO']
    # create_eventclasses = ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']
    # create_eventclasses = ['OIBEO', 'IBEO', 'OIBE']

    submit_eventclasses = ['OIBEO']
    # submit_eventclasses = ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']
    # submit_eventclasses = ['OIBEO', 'IBEO', 'OIBE']
    # submit_eventclasses = ['OIB', 'BEO']

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

if __name__ == "__main__":
    main()
