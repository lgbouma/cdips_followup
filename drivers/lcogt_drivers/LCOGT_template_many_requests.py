"""
"I am in location X. Is there any CDIPS object I can observe tonight?"
"""
import numpy as np
from astropy.time import Time
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.LCOGT_dedicated_requests import (
    get_dedicated_request,
    given_dedicated_requests_validate_submit
)

if __name__ == "__main__":

    ####################
    overwrite = 1
    validate = 0
    submit = 0

    savstr = '20191228_bakos_special' # eg, 20191207_TOI1098_request_2m_tc_secondary. "ephemupdate" if it is one. (this cancels pending observations)

    source_ids = [
        4844691297067063424,
        5519619186857962112,
        5525188767305211904,
        5251470948229949568,
        5290721997195236480,
        5576476552334683520,
        5557593814516968960,
        5579734916388215808,
        5325454783543157760,
        5838450865699668736,
        2919143383943171200,
        3340674976430098688,
        5290781443841554432,
        5510676828723793920,
        5516140233292943872,
        5561614350584396800,
        5596735638203997824
    ]
    source_ids = np.array(source_ids).astype(str)

    max_search_time = Time('2019-12-30 23:59:00')
    filtermode = 'ip'
    telescope_class = '1m0'

    # submit_eventclasses = ['OIBEO' 'IBEO', 'OIBE', 'OIB', 'BEO']
    # submit_eventclasses = ['OIBEO']
    # submit_eventclasses = ['OIB', 'BEO']
    submit_eventclasses = [] # ['OIBEO', 'IBEO', 'OIBE']

    manual_ephemeris = False
    max_duration_error = 20
    create_eventclasses = ['OIB', 'BEO'] # 'OIBEO', 'IBEO', 'OIBE', 
    ####################

    for source_id in source_ids:

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
