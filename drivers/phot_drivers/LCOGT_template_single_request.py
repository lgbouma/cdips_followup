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
    savstr = '20250217_toi1224_25b' # eg, 20191207_TOI1098_request_2m_tc_secondary. "ephemupdate" if it is one. (this cancels pending observations)
    overwrite = 1
    validate = 0
    submit = 0

    tic_id = '299798795' # '120105470'
    source_id = None # '6113920619134019456' # can use instead of TIC

    filtermode = 'ip'# 'zs', 'gp', 'ip'
    #telescope_class = '1m0' # '1m0', '2m0', 'special'
    telescope_class = 'special' # '1m0', '2m0', 'special'
    ipp_value = 1 # usually 1
    #max_search_time = Time('2022-12-31 23:59:00')
    max_search_time = Time('2027-01-31 23:59:00')

    verify_ephemeris_uncertainty = 1 # require t_tra uncertainty < 2 hours
    inflate_duration = 0 # if t_tra uncertainty > 1 hour, inflate tdur by +/- 45 minutes per side

    transit_type = 'totals' # 'totals', 'partials', 'all', 'fulltotals' -- see above
    max_n_events = 99 # else None. n_events is per eventclass.

    raise_error = False # raise an error if max_duration_error flag raised.
    max_duration_error = 30 # the submitted LCOGT request must match requested durn to within this difference [minutes]
    #sites = ['Palomar'] # Default None for LCOGT. Could do e.g., 'special' and ['Keck Observatory']
    sites = ['Las Campanas Observatory'] # Default None for LCOGT. Could do e.g., 'special' and ['Keck Observatory']
    #sites = ['Whipple Observatory'] # Default None for LCOGT. Could do e.g., 'special' and ['Keck Observatory']
    #sites = ['Keck Observatory'] # Default None for LCOGT. Could do e.g., 'special' and ['Keck Observatory']
    #sites = ['Cerro Paranal'] # Default None for LCOGT. Could do e.g., 'special' and ['Keck Observatory']

    force_acceptability = 50 # None or int.
    semesterstr = '25B'

    # CHANGE ABOVE
    ##########################################

    max_airmass_sched = 2.5
    manual_ephemeris = False
    manual_ephemeris = True # FIXME
    create_eventclasses = TRANSITTYPEDICT[transit_type]
    submit_eventclasses = TRANSITTYPEDICT[transit_type]

    if source_id is None:
        assert isinstance(tic_id, str)
        source_id = tic_to_gaiadr2(tic_id)

    if manual_ephemeris:
        period = 17.945466
        period_unc = 0.000012
        epoch = 2458329.859704
        #epoch = 2460260.25 # half an orbit later
        #epoch = 2460286.71 # the new dip
        epoch_unc = 0.000828
        duration = 3

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
            msg = (f'ERR! Got negative ephem unc of '+
                   f'{delta_t_tra_today*24:.1f} hr. '+
                   f'Need to give a believable ephem unc..')
            raise ValueError(msg)
        if delta_t_tra_today*24 > 2:
            msg = f'ERR! Got ephem unc of {delta_t_tra_today*24:.3f} hr. This is too high.'
            raise ValueError(msg)
        if delta_t_tra_today*24 > 1:
            msg = f'WRN! Got ephem unc of {delta_t_tra_today*24:.3f} hr. This is risky.'
            print(msg)
        else:
            msg = f'INFO! Got ephem unc of {delta_t_tra_today*24:.3f} hr. This is fine.'
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
        semesterstr=semesterstr,
        overwrite=overwrite, max_search_time=max_search_time,
        filtermode=filtermode, telescope_class=telescope_class,
        ipp_value=ipp_value, sites=sites,
        force_acceptability=force_acceptability,
        max_airmass_sched=max_airmass_sched
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
