"""
Given targets, their ephemerides, and positions, create requests to LCOGT to
get imaging follow-up with the 1m and 2m.
"""

###########
# imports #
###########
import requests, socket, os, pickle
import numpy as np, pandas as pd
from numpy import array as nparr

from get_event_observability import \
        get_event_observability, print_event_observability

import datetime as dt
from astropy.time import Time
from astropy.coordinates import get_body, get_sun, get_moon, SkyCoord
import astropy.units as u

from astroplan import (FixedTarget, Observer, EclipsingSystem,
                       PrimaryEclipseConstraint, is_event_observable,
                       AtNightConstraint, AltitudeConstraint,
                       LocalTimeConstraint, MoonSeparationConstraint,
                       moon)

from astroquery.gaia import Gaia

if socket.gethostname() == 'brik':
    api_file = '/home/luke/.lcogt_api_token'
elif 'astro' in socket.gethostname():
    api_file = '/Users/luke/.lcogt_api_token'
else:
    raise NotImplementedError('where to get API file?')

with open(api_file, 'r') as f:
    l = f.readlines()
token = str(l[0].replace('\n',''))

#############
# functions #
#############

def _given_Gmag_get_exptime_defocus(Gmag, telescope_class):
    """
    this table was reverse-engineered by looking at the exptime and defocuses
    used by Dan Bayliss in the HATS requests for the same proposal.

    it imposes a minimum time 10 second exposures, or a bright limit of Gaia G
    mag of G=9.

    (any brighter and we might encounter smear issues -- need to test).
    """
    if telescope_class == '1m0':
        df = pd.read_csv('../data/LCOGT_reverse_eng_exptime.csv')
    elif telescope_class == '2m0':
        df = pd.read_csv('../data/LCOGT_2m_guess_exptime.csv')

    if Gmag > df.G.max():
        raise AssertionError('target too faint')
    elif Gmag < df.G.min():
        raise AssertionError('target too bright')

    temp = df.iloc[(df['G']-Gmag).abs().argsort()[0]]

    exptime = temp['exptime']
    defocus = temp['defocus']

    return exptime, defocus


def make_request_group(targetname, ra, dec, pmra, pmdec, Gmag, starttime,
                       endtime, eventclass='OIBEO', max_airmass=2.5,
                       min_lunar_distance=20, filtermode="ip",
                       telescope_class="1m0"):

    try:
        exptime, defocus = _given_Gmag_get_exptime_defocus(
            Gmag, telescope_class
        )
    except AssertionError as e:
        print(e)
        return -1

    API_TOKEN = token  # API token obtained from https://observe.lco.global/accounts/profile/
    PROPOSAL_ID = 'NOAO2019B-013'  # Proposal IDs may be found here: https://observe.lco.global/proposals/

    # starttime e.g., '2019-05-02 00:00:00'
    _starttime = starttime.iso[0:19]
    _endtime = endtime.iso[0:19]

    read_time_per_exposure = 30*u.second # from Bayliss' completed runs

    expcount = np.floor(
        (endtime-starttime).to(u.hr)
        /
        (exptime*u.second + read_time_per_exposure).to(u.hr)
    )

    obsdate = (
          str(_starttime).split('-')[0]      # year
        + str(_starttime).split('-')[1]      # month
        + str(_starttime).split('-')[2][:2]  # date
    )

    requestname = (
        '{targetname:s}_{eventclass:s}_{telescope_class:s}_{obsdate:s}_{filtermode:s}_{exptime:d}_{defocus:.1f}'.
        format(targetname=targetname, eventclass=eventclass,
               telescope_class=telescope_class, obsdate=obsdate,
               filtermode=filtermode, exptime=exptime, defocus=defocus)
    )

    # The target of the observation
    target = {
        'name': str(targetname),
        'type': 'ICRS',
        'ra': float(ra), # decimal deg
        'dec': float(dec),
        'proper_motion_ra': float(pmra),
        'proper_motion_dec': float(pmdec),
        'epoch': 2015.5
    }

    # Constraints used for scheduling the observation
    constraints = {
        'max_airmass': max_airmass,
        'min_lunar_distance': min_lunar_distance
    }

    # The configurations for this request. In this example we are taking 2 exposures with
    # different filters and exposure times. The fields acquisition_config and guiding_config 
    # are required fields in a configuration that are ultimately filled in with defaults 
    # if the submitted values are empty.
    if telescope_class == '1m0':
        instrument_type = '1M0-SCICAM-SINISTRO'
        mode = 'full_frame'
        bin_x, bin_y = 1, 1
    elif telescope_class == '2m0':
        instrument_type = '2M0-SCICAM-SPECTRAL'
        mode = 'default'
        bin_x, bin_y = 2, 2

    configurations = [
        {
            'type': 'EXPOSE',
            'instrument_type': instrument_type,
            'target': target,
            'constraints': constraints,
            "instrument_configs": [
                {
                    "optical_elements": {
                        "filter": filtermode
                    },
                    "mode": mode,
                    "exposure_time": int(exptime),
                    "exposure_count": int(expcount),
                    "bin_x": bin_x,
                    "bin_y": bin_y,
                    "rotator_mode": "",
                    "extra_params": {
                        "defocus": float(defocus)
                    }
                }
            ],
            'acquisition_config': {},
            'guiding_config': {
                "optional": False,
                "mode": "ON"
            },
            "priority": 1
        }
    ]

    # The time windows during which this request should be considered for observing. In this example
    # we only provide one. These times are in UTC.
    windows = [{
        'start': _starttime, # e.g., '2019-05-02 00:00:00',
        'end': _endtime # '2019-05-30 00:00:00'
    }]

    # The telescope class that should be used for this observation
    location = {
        'telescope_class': telescope_class
    }

    # The full RequestGroup, with additional meta-data
    requestgroup = {
        'name': requestname,  # The title
        'proposal': PROPOSAL_ID,
        'ipp_value': 1.05,
        'operator': 'SINGLE',
        'observation_type': 'NORMAL',
        'requests': [{
            'configurations': configurations,
            'windows': windows,
            'location': location,
        }]
    }

    return requestgroup


def get_requests_given_ephem(
    savstr, targetname, ra, dec, pmra, pmdec, Gmag, period, period_unc, epoch,
    epoch_unc, depth, depth_unc, duration, duration_unc,
    min_search_time=Time(dt.datetime.today().isoformat()),
    max_search_time=Time('2019-11-30 23:59:00'),
    max_airmass_sched=1.8,
    max_airmass_submit=2.5,
    min_lunar_distance=20, oot_duration=45*u.minute,
    eventclass='OIBEO',
    sites=['Cerro Tololo', 'Siding Spring Observatory', 'SAAO'],
    schedule_oot_duration=60*u.minute):
    """
    Given an ephemeris, and the basic details of a target, generate LCOGT
    requests for any available transits at the given sites, between
    min_search_time and max_search_time.

    Allowed sites include Siding Spring Observator and CTIO.

    Args:

        get_oibeo: gets the request for which "OIBEO" transits visible (given
        oot_duration).

        get_ibe: just "IBE" transit required (easier to schedule).

        schedule_oot_duration: used for the LCOGT-side REQUESTS, rather than
        the astroplan scheduling. Can be longer than oot_duration (used for
        astroplan scheduling) in order to try and get a bit more data.

    Note:
        LCO Semester A is Dec 1 thru May 31.
        LCO Semester B is June 1 thru Nov 30.
    """
    if eventclass not in [
        'OIBEO', 'OIBE', 'IBEO', 'IBE', 'BEO', 'OIB'
    ]:
        raise AssertionError

    outdir = "../results/LCOGT_19B20A_observability/{}".format(savstr)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outdir = os.path.join(outdir, '{}'.format(targetname))
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if max_airmass_sched>3:
        raise NotImplementedError('approx breaks')

    groups, sel_sites = [], []

    for site in sites:

        _site = Observer.at_site(site)

        (event_ind, oibeo,
         ing_tmid_egr, target_window,
         moon_separation, moon_illumination
        ) = (
            get_event_observability(
                eventclass,
                _site, ra*u.deg, dec*u.deg, targetname, epoch, period*u.day,
                duration*u.hour, n_transits=100,
                obs_start_time=min_search_time,
                oot_duration=oot_duration,
                minokmoonsep=min_lunar_distance*u.deg,
                max_airmass=max_airmass_sched
            )
        )

        lt_maxtime = ing_tmid_egr < np.array(max_search_time)

        sel = event_ind.flatten() & np.all(lt_maxtime, axis=1)

        sel_times = target_window[sel, :]

        if len(sel_times)>=1:

            ##############################################
            # save observability output as text file too #
            print_sel = (
                event_ind.flatten()
                &
                np.all(lt_maxtime, axis=1)
            )
            print_event_observability(eventclass,
                                      event_ind[:,print_sel],
                                      oibeo[:,print_sel],
                                      ing_tmid_egr[print_sel,:], _site,
                                      ra*u.deg, dec*u.deg, targetname, epoch,
                                      period*u.day, duration*u.hour,
                                      max_airmass_sched, oot_duration,
                                      moon_separation[print_sel],
                                      moon_illumination[print_sel],
                                      minokmoonsep=min_lunar_distance*u.deg,
                                      outdir=outdir)
            ##############################################

            if '_2m_' in savstr:
                telescope_class = "2m0"
            else:
                telescope_class = "1m0"

            for sel_time in sel_times:

                assert len(sel_time) == 2

                starttime = sel_time[0]
                endtime = sel_time[1]

                g = make_request_group(targetname, ra, dec, pmra,
                                       pmdec, Gmag, starttime, endtime,
                                       eventclass=eventclass, filtermode="ip",
                                       telescope_class=telescope_class,
                                       max_airmass=max_airmass_submit,
                                       min_lunar_distance=min_lunar_distance)

                if g == -1:
                    continue
                else:
                    groups.append(g)

    return groups


def get_all_requests_19B(savstr, eventclass):

    df = get_targets(savstr, verbose=False)

    results = []

    for ix, r in df.iterrows():

        #
        # get gaia positions and PMs (the coordinates read in are slightly off)
        #
        jobstr = (
            "select top 1 g.ra, g.dec, g.pmra, g.pmdec from "
            "gaiadr2.gaia_source as g where g.source_id = {:d}".
            format(np.int64(r['source_id']))
        )

        job = Gaia.launch_job(jobstr)
        gaia_r = job.get_results()

        if len(gaia_r) != 1:
            raise AssertionError('gaia match failed')

        ra = float(gaia_r['ra'])
        dec = float(gaia_r['dec'])
        pmra = float(gaia_r['pmra'])
        pmdec = float(gaia_r['pmdec'])

        #
        # shift by 42 arcseconds away from the center, in order to avoid CCD
        # amplifier lines.
        #
        c = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')

        shift_by = 42*u.arcsec # Bayliss shifted by ~30 arcsec. might as well further.
        shift_dir = 45*u.deg   # as long as it's some mix of "up" and "left"

        use_coord = c.directional_offset_by(shift_dir, shift_by)
        ra = use_coord.ra.value
        dec = use_coord.dec.value

        #
        #
        #
        if '_2m_' in savstr:
            sites = ['Siding Spring Observatory', 'Haleakala Observatories']
        else:
            # assume 1m
            sites = ['Cerro Tololo', 'Siding Spring Observatory', 'SAAO']

        this = get_requests_given_ephem(savstr, r['toi_or_ticid'],
                                        ra, dec, pmra, pmdec,
                                        r['phot_g_mean_mag'], r['period'],
                                        r['period_unc'], r['epoch'],
                                        r['epoch_unc'], r['depth'],
                                        r['depth_unc'], r['duration'],
                                        r['duration_unc'], sites=sites,
                                        eventclass=eventclass)

        results.append(this)

    return results


def main(savstr=None, overwrite=None, eventclass=None):
    """
    savstr:

        'all_requests_19B_easyones': original request, G=9-15.4, depth>500
        'request_19B_2m_faint': two faint boyos for the 2m. (one didnt work)
        'request_19B_2m_faint_v2': getting the third faint (PMS) one
        'request_TIC29786532_19B': on the 1m, schedule the one that didnt work

    eventclass:
        any of "OIBEO", "IBEO", "BEO", etc.
    """

    assert isinstance(savstr, str)
    assert isinstance(overwrite, int)

    pkl_savpath = (
        '../results/LCOGT_19B20A_observability/{}.pkl'.format(savstr)
    )
    mult_savpath = (
        '../results/LCOGT_19B20A_observability/{}_summary.csv'.format(savstr)
    )

    if not overwrite and os.path.exists(pkl_savpath):
        with open(pkl_savpath, 'rb') as f:
            r = pickle.load(f)
    else:
        r = get_all_requests_19B(savstr, eventclass)
        with open(pkl_savpath, 'wb') as f:
            pickle.dump(r, f, pickle.HIGHEST_PROTOCOL)
            print('saved {:s}'.format(pkl_savpath))

    df = get_targets(savstr, verbose=True)

    names = [_r['toi_or_ticid'] for ix, _r in df.iterrows()]
    mult = [len(_r) for _r in r]
    durns = [_r['duration'] for ix, _r in df.iterrows()]
    durns_sched = [_r['duration']+2 for ix, _r in df.iterrows()]
    mult_df = pd.DataFrame({
        'name':names,
        'n_requests':mult,
        'duration':durns,
        'sched_duration':durns_sched,
        'n*sched':np.array(mult)*np.array(durns_sched)
    })
    mult_df.to_csv(mult_savpath, index=False)
    print('made {}'.format(mult_savpath))

    return r, mult_df


def get_targets(savstr, verbose=True):

    df = pd.read_csv('../data/20190912_19B20A_LCOGT_1m_2m.csv')

    if savstr in ['all_requests_19B_easyones']:
        sel = (
            (df['phot_g_mean_mag'] < 15.4)
            &
            (df['phot_g_mean_mag'] > 9)
            &
            (df['depth'] > 500) # 500 ppm = 0.05% = 0.5 mmag
        )

    elif savstr == 'request_19B_2m_faint':
        ids = ['TIC29786532.01', 'TIC53682439.01'] # faint
    elif savstr == 'request_19B_2m_faint_v2':
        ids = ['TIC200516835.01'] # faint
    elif savstr == 'request_TIC29786532_19B':
        ids = ['TIC29786532.01'] # faint
    elif savstr in ['toppartials_19B_IBEO', 'toppartials_19B_OIBE',
                    'toppartials_19B_IBE', 'toppartials_19B_BEO',
                    'toppartials_19B_OIB']:
        ids = ['TIC238611475.01',
               'TIC125192758.01',
               'TIC154671430.01',
               'TIC110718787.01', # these 4 have no OIBEO transits
               'TIC308538095.01' # only 1 OIBEO
              ]
    else:
        raise NotImplementedError

    if savstr in ['request_19B_2m_faint', 'request_19B_2m_faint_v2',
                  'request_TIC29786532_19B', 'toppartials_19B_IBEO',
                  'toppartials_19B_OIBE', 'toppartials_19B_IBE',
                  'toppartials_19B_BEO', 'toppartials_19B_OIB']:
        sel = (
            df['toi_or_ticid'].isin(ids)
        )

    newdf = df[sel]

    if verbose:
        print(42*'-')
        print('WRN: REQUEST WAS {}'.format(savstr))
        print('WRN: DROPPING THE FOLLOWING TARGETS B/C OUTSIDE DESIRED REQUEST')
        print(df[~sel][['source_id', 'toi_or_ticid']])
        print(42*'-')

    return newdf



if __name__ == "__main__":
    #NOTE : "_2m_" must be in savstr to request 2m settings

    eventclass = 'IBEO'
    savstr = 'toppartials_19B_{}'.format(eventclass)

    ##########
    # eventclass = 'OIBEO'
    # savstr = 'request_TIC29786532_19B'
    # savstr = 'request_19B_2m_faint_v2'
    # savstr = 'request_19B_2m_faint'
    # savstr = 'all_requests_19B_easyones'
    ##########

    overwrite = 1

    r, mult_df = main(savstr=savstr, overwrite=overwrite, eventclass=eventclass)
