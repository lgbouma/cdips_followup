from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.LCOGT_dedicated_requests import (
    get_dedicated_request,
    given_dedicated_requests_validate_submit
)

if __name__ == "__main__":

    ####################
    overwrite = 0
    validate = 1
    submit = 0
    max_duration_error = 20
    manual_ephemeris = False
    savstr = '20191207_TOI837_request_1m' # eg, 20191207_TOI1098_request_2m_tc_secondary. "ephemupdate" if it is one.
    source_id = '5251470948229949568'
    ####################

    if manual_ephemeris:
        period = 42
        epoch = 2458660.00000
        duration = 2.00000
    else:
        # get ephemeris from ephermides.csv
        d = query_ephemeris(source_id=source_id)
        period, epoch, duration = d['period'], d['epoch'], d['duration']

    eventclasses = ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO']

    requests = get_dedicated_request(savstr, source_id, period, epoch,
                                     duration, eventclasses,
                                     overwrite=overwrite)

    given_dedicated_requests_validate_submit(
        requests, validate=validate, submit=submit,
        max_duration_error=max_duration_error)
