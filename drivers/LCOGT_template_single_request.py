from astropy.time import Time
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.LCOGT_dedicated_requests import (
    get_dedicated_request,
    given_dedicated_requests_validate_submit
)

if __name__ == "__main__":

    ####################
    overwrite = 0
    validate = 1
    submit = 1
    max_duration_error = 20
    savstr = '20191209_1m_requests' # eg, 20191207_TOI1098_request_2m_tc_secondary. "ephemupdate" if it is one.
    source_id = '5596735638203997824'
    manual_ephemeris = False
    max_search_time = Time('2020-01-15 23:59:00')

    eventclasses = ['OIBEO'] #['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']
    ####################

    if manual_ephemeris:
        period = 42
        epoch = 2458660.00000
        duration = 2.00000

    else:
        # get ephemeris from ephermides.csv
        d = query_ephemeris(source_id=source_id)
        period, epoch, duration = d['period'], d['epoch'], d['duration']

    # "requests" is a list of lists. Higher level is each eventclass. Level
    # below is each event, in that eventclass.
    requests = get_dedicated_request(savstr, source_id, period, epoch,
                                     duration, eventclasses,
                                     overwrite=overwrite,
                                     max_search_time=max_search_time)

    given_dedicated_requests_validate_submit(
        requests, validate=validate, submit=submit,
        max_duration_error=max_duration_error)
