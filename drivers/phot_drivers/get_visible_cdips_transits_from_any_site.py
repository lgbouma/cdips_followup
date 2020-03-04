"""
"I am in location X.
Is there any CDIPS object I can observe between time Y and Z?"

The solution here hacks the LCOGT drivers to answer this question.
"""

import numpy as np
from astropy.time import Time
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.LCOGT_dedicated_requests import (
    get_dedicated_request,
    given_dedicated_requests_validate_submit
)
from cdips_followup.manage_candidates import get_live_candidates

if __name__ == "__main__":

    ####################
    # vvv change below
    savstr = '20200304_deLeon_special'
    condition = 'SP0SP1'

    overwrite = 0
    min_search_time = Time('2020-03-22 11:59:00')
    max_search_time = Time('2020-04-01 23:59:00')

    sites=['SAAO']
    create_eventclasses = ['OIBEO']
    # ^^^ change above
    ####################

    # sites=['Cerro Tololo', 'Siding Spring Observatory', 'SAAO',
    #        'McDonald Observatory']
    filtermode = 'ip'
    telescope_class = 'special'

    submit_eventclasses = []
    max_duration_error = 20
    validate = 0
    submit = 0
    ####################

    source_ids = get_live_candidates(condition=condition)
    source_ids = np.array(source_ids).astype(str)

    for source_id in source_ids:

        # get ephemeris from ephemerides.csv
        try:
            d = query_ephemeris(source_id=source_id)
        except Exception as e:
            print('WRN! Failed to get ephemeris for {}'.format(source_id))
            continue
        period, epoch, duration = d['period'], d['epoch'], d['duration']

        # "requests" is a list of lists. Higher level is each eventclass. Level
        # below is each event, in that eventclass.
        requests = get_dedicated_request(savstr, source_id, period, epoch,
                                         duration, create_eventclasses,
                                         overwrite=overwrite,
                                         min_search_time=min_search_time,
                                         max_search_time=max_search_time,
                                         filtermode=filtermode,
                                         telescope_class=telescope_class,
                                         sites=sites)

        given_dedicated_requests_validate_submit(
            requests, submit_eventclasses, validate=validate, submit=submit,
            max_duration_error=max_duration_error)
