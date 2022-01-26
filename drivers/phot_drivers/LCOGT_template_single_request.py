"""
Given a source_id, make LCOGT photometry followup requests, and optionally
submit them to the LCOGT API.
"""
import numpy as np
from astropy.time import Time
from cdips_followup.manage_ephemerides import (
    query_ephemeris, get_ephemeris_uncertainty
)
from cdips_followup.LCOGT_dedicated_requests import (
    get_dedicated_request,
    given_dedicated_requests_validate_submit
)
from astrobase.services.identifiers import tic_to_gaiadr2

TRANSITTYPEDICT = {
    'all': ['OIBEO', 'IBEO', 'OIBE', 'OIB', 'BEO'],
    'partials': ['OIB', 'BEO'],
    'totals': ['OIBEO', 'IBEO', 'OIBE'],
    'fulltotals': ['OIBEO']
}

def main():

    ##########################################
    # CHANGE BELOW
    savstr = '20211211_v1096tau' # eg, 20191207_TOI1098_request_2m_tc_secondary. "ephemupdate" if it is one. (this cancels pending observations)
    overwrite = 1
    validate = 1
    submit = 1

    tic_id = '56551765' # '120105470'
    source_id = None # '6113920619134019456' # can use instead of TIC

    filtermode = 'rp'# 'zs', 'gp', 'ip'
    telescope_class = '1m0' # '1m0', '2m0', 'special'
    #telescope_class = 'special' # '1m0', '2m0', 'special'
    ipp_value = 1 # usually 1
    #max_search_time = Time('2022-12-31 23:59:00')
    max_search_time = Time('2022-01-31 23:59:00')

    verify_ephemeris_uncertainty = 1 # require t_tra uncertainty < 2 hours
    inflate_duration = 1 # if t_tra uncertainty > 1 hour, inflate transit duration by +/- 45 minutes per side

    transit_type = 'totals' # see above
    max_n_events = 99 # else None. n_events is per eventclass.

    raise_error = False # raise an error if max_duration_error flag raised.
    max_duration_error = 30 # the submitted LCOGT request must match requested durn to within this difference [minutes]

    sites = None #['Keck Observatory'] # Default None for LCOGT. Could do e.g., 'special' and ['Keck Observatory']
    #sites = ['Keck Observatory'] # Default None for LCOGT. Could do e.g., 'special' and ['Keck Observatory']

    force_acceptability = 50 # None or int.

    # CHANGE ABOVE
    ##########################################

    manual_ephemeris = False
    create_eventclasses = TRANSITTYPEDICT[transit_type]
    submit_eventclasses = TRANSITTYPEDICT[transit_type]

    if source_id is None:
        assert isinstance(tic_id, str)
        source_id = tic_to_gaiadr2(tic_id)

    if manual_ephemeris:
        period = 42
        epoch = 2458660.00000
        duration = 2.00000

    else:
        # get ephemeris from ephemerides.csv
        d = query_ephemeris(source_id=source_id)
        period, epoch, duration = (
            d['period'], d['epoch'], d['duration']
        )
        period_unc, epoch_unc, duration_unc = (
            d['period_unc'], d['epoch_unc'], d['duration_unc']
        )

    if verify_ephemeris_uncertainty:
        delta_t_tra_today = (
            get_ephemeris_uncertainty(epoch, epoch_unc, period, period_unc, epoch_obs='today')
        )
        if delta_t_tra_today*24 < 0:
            msg = f'ERR! Got negative ephem unc of {delta_t_tra_today*24:.1f} hr. Need to give a believable ephem unc..'
            raise ValueError(msg)
        if delta_t_tra_today*24 > 2:
            msg = f'ERR! Got ephem unc of {delta_t_tra_today*24:.1f} hr. This is too high.'
            raise ValueError(msg)
        if delta_t_tra_today*24 > 1:
            msg = f'WRN! Got ephem unc of {delta_t_tra_today*24:.1f} hr. This is risky.'
            print(msg)

    if inflate_duration:
        assert verify_ephemeris_uncertainty
        if delta_t_tra_today*24 > 1:
            msg = f'... inflating transit duration for scheduling pursposes by 1.5 hours.'
            print(msg)
            duration += 1.5 # add

    # "requests" is a list of lists. Higher level is each eventclass. Level
    # below is each event, in that eventclass.
    requests = get_dedicated_request(
        savstr, source_id, period, epoch, duration, create_eventclasses,
        overwrite=overwrite, max_search_time=max_search_time,
        filtermode=filtermode, telescope_class=telescope_class,
        ipp_value=ipp_value, sites=sites,
        force_acceptability=force_acceptability
    )

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

        if len(_requests) == 0:
            print('WRN!: got no times')
            return

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
        max_duration_error=max_duration_error, raise_error=raise_error
    )


if __name__ == "__main__":
    main()
