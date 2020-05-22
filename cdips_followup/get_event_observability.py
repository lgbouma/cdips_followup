import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
from glob import glob

from astropy.time import Time
from astropy.coordinates import get_body, get_sun, get_moon, SkyCoord
import astropy.units as u

from astroplan import (FixedTarget, Observer, EclipsingSystem,
                       PrimaryEclipseConstraint, is_event_observable,
                       AtNightConstraint, AltitudeConstraint,
                       LocalTimeConstraint, MoonSeparationConstraint,
                       AirmassConstraint, moon)

import datetime as dt


def get_event_observability(
    eventclass,
    site, ra, dec, name, t_mid_0, period, duration, n_transits=100,
    obs_start_time=Time(dt.datetime.today().isoformat()),
    min_altitude = None,
    oot_duration = 30*u.minute,
    minokmoonsep = 30*u.deg,
    max_airmass = None,
    twilight_limit = 'nautical'):
    """
    note: barycentric corrections not yet implemented. (could do this myself!)
    -> 16 minutes of imprecision is baked into this observability calculator!

    args:

        eventclass: e.g., "OIBE". Function does NOT return longer events.

        site (astroplan.observer.Observer)

        ra, dec (units u.deg), e.g.:
            ra=101.28715533*u.deg, dec=16.71611586*u.deg,
        or can also accept
            ra="17 56 35.51", dec="-29 32 21.5"

        name (str), e.g., "Sirius"

        t_mid_0 (float): in BJD_TDB, preferably (but see note above).

        period (astropy quantity, units time)

        duration (astropy quantity, units time)

        n_transits (int): number of transits forward extrapolated to

        obs_start_time (astropy.Time object): when to start calculation from

        min_altitude (astropy quantity, units deg): 20 degrees is the more
        relevant constraint.

        max_airmass: e.g., 2.5. One of max_airmass or min_altitude is required.

        oot_duration (astropy quantity, units time): with which to brack
        transit observations, to get an OOT baseline.

        twilight_limit: 'astronomical', 'nautical', 'civil' for -18, -12, -6
        deg.
    """
    if eventclass not in [
        'OIBEO', 'OIBE', 'IBEO', 'IBE', 'BEO', 'OIB', 'OI', 'EO'
    ]:
        raise AssertionError

    if (isinstance(ra, u.quantity.Quantity) and
        isinstance(dec, u.quantity.Quantity)
    ):
        target_coord = SkyCoord(ra=ra, dec=dec)
    elif (isinstance(ra, str) and
          isinstance(dec, str)
    ):
        target_coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
    else:
        raise NotImplementedError

    if (
        not isinstance(max_airmass, float)
        or isinstance(min_altitude, u.quantity.Quantity)
    ):
        raise NotImplementedError

    target = FixedTarget(coord=target_coord, name=name)

    primary_eclipse_time = Time(t_mid_0, format='jd')

    system = EclipsingSystem(primary_eclipse_time=primary_eclipse_time,
                             orbital_period=period, duration=duration,
                             name=name)

    midtransit_times = system.next_primary_eclipse_time(
        obs_start_time, n_eclipses=n_transits)

    # for the time being, omit any local time constraints.
    if twilight_limit == 'astronomical':
        twilight_constraint = AtNightConstraint.twilight_astronomical()
    elif twilight_limit == 'nautical':
        twilight_constraint = AtNightConstraint.twilight_nautical()
    else:
        raise NotImplementedError('civil twilight is janky.')

    constraints = [twilight_constraint,
                   AltitudeConstraint(min=min_altitude),
                   AirmassConstraint(max=max_airmass),
                   MoonSeparationConstraint(min=minokmoonsep)]

    # tabulate ingress and egress times.
    ing_egr = system.next_primary_ingress_egress_time(
        obs_start_time, n_eclipses=n_transits
    )

    oibeo_window = np.concatenate(
        (np.array(ing_egr[:,0] - oot_duration)[:,None],
         np.array(ing_egr[:,1] + oot_duration)[:,None]),
        axis=1)
    oibe_window = np.concatenate(
        (np.array(ing_egr[:,0] - oot_duration)[:,None],
         np.array(ing_egr[:,1])[:,None]),
        axis=1)
    ibeo_window = np.concatenate(
        (np.array(ing_egr[:,0])[:,None],
         np.array(ing_egr[:,1] + oot_duration)[:,None]),
        axis=1)
    oib_window = np.concatenate(
        (np.array(ing_egr[:,0] - oot_duration)[:,None],
         np.array(midtransit_times)[:,None]),
        axis=1)
    beo_window = np.concatenate(
        (np.array(midtransit_times)[:,None],
         np.array(ing_egr[:,1] + oot_duration)[:,None]),
        axis=1)
    ibe_window = ing_egr
    oi_window = np.concatenate(
        (np.array(ing_egr[:,0] - oot_duration)[:,None],
        np.array(ing_egr[:,0])[:,None]),
        axis=1)
    eo_window = np.concatenate(
        (np.array(ing_egr[:,1])[:,None],
        np.array(ing_egr[:,1] + oot_duration)[:,None]),
        axis=1)

    keys = ['oibeo','oibe','ibeo','oib','beo','ibe','oi','eo']
    windows = [oibeo_window, oibe_window, ibeo_window,
               oib_window, beo_window, ibe_window, oi_window, eo_window]
    is_obs_dict = {}
    for key, window in zip(keys, windows):
        is_obs_dict[key] = np.array(
            is_event_observable(constraints, site, target,
                                times_ingress_egress=window)
        ).flatten()

    is_obs_df = pd.DataFrame(is_obs_dict)

    is_obs_df['ing'] = ing_egr[:,0]
    is_obs_df['egr'] = ing_egr[:,1]
    is_obs_df['isoing'] = Time(ing_egr[:,0], format='iso')
    is_obs_df['isoegr'] = Time(ing_egr[:,1], format='iso')

    # this function returns the observable events that are LONGEST. e.g.,
    # during an OIBEO transit you COULD observe just OIB, but why would you?

    if eventclass == 'OIBEO':
        event_ind = np.array(is_obs_df[eventclass.lower()])[None,:]
    elif eventclass in ['IBEO', 'OIBE']:
        event_ind = np.array(
            is_obs_df[eventclass.lower()] & ~is_obs_df['oibeo']
        )[None,:]
    elif eventclass in ['IBE', 'OIB', 'BEO']:
        event_ind = np.array(
            is_obs_df[eventclass.lower()]
            & ~is_obs_df['oibeo']
            & ~is_obs_df['oibe']
            & ~is_obs_df['ibeo']
        )[None,:]
    elif eventclass in ['OI', 'EO']:
        event_ind = np.array(
            is_obs_df[eventclass.lower()]
            & ~is_obs_df['oibeo']
            & ~is_obs_df['oibe']
            & ~is_obs_df['ibeo']
            & ~is_obs_df['oib']
            & ~is_obs_df['ibe']
            & ~is_obs_df['beo']
        )[None,:]

    # get moon separation over each transit. take minimum moon sep at
    # ing/tmid/egr as the moon separation.
    moon_tmid = get_moon(midtransit_times, location=site.location)
    moon_separation_tmid = moon_tmid.separation(target_coord)

    moon_ing = get_moon(ing_egr[:,0], location=site.location)
    moon_separation_ing = moon_ing.separation(target_coord)

    moon_egr = get_moon(ing_egr[:,1], location=site.location)
    moon_separation_egr = moon_egr.separation(target_coord)

    moon_separation = np.round(np.array(
        [moon_separation_tmid, moon_separation_ing,
         moon_separation_egr]).min(axis=0),0).astype(int)

    moon_illumination = np.round(
        100*moon.moon_illumination(midtransit_times),0).astype(int)

    # completely observable transits (OOT, ingress, bottom, egress, OOT)
    oibeo = is_event_observable(constraints, site, target,
                                times_ingress_egress=oibeo_window)

    ing_tmid_egr = np.concatenate(
        (np.array(ing_egr[:,0])[:,None],
         np.array(midtransit_times)[:,None],
         np.array(ing_egr[:,1])[:,None]),
        axis=1)

    target_window = np.array(windows)[
        int(np.argwhere(np.array(keys)==eventclass.lower())), :, :
    ]

    return (
        event_ind, oibeo, ing_tmid_egr, target_window,
        moon_separation, moon_illumination
    )


def print_event_observability(eventclass, event_ind, oibeo, ing_tmid_egr, site,
                              ra, dec, name, t_mid_0, period, duration,
                              max_airmass, oot_duration, moon_separation,
                              moon_illumination, minokmoonsep=30*u.deg,
                              outdir="../results/"):
    """
    write output to ../results/{name}_{site}.txt
    """

    # header information
    hdr = (
"""{} transit prediction

site = {}.
lon {:.2f}, lat {:.2f}, alt {:.1f}

RA = {}, DEC = {}.

NOTE: This calculator does not yet have barycenter correction implemented! This
means that all local times of astronomical events are at worst 16 minutes off.

Given:
        t0_BJD_TDB = {}
        period = {}
        duration = {}

Constrain to observe with:
        airmass < {}
        and
        night-time (nautical twilight)
        and
        moon separation > {}

Defined "out of transit" to be {} on either side of ingress/egress.

==========================================
OBSERVABLE IN {} (ALL TIMES IN UTC)

epoch     t_ing                 t_mid                 t_egr                 sep  illum
--------------------------------------------------------------------------------------
""".format(
    name, site.name,
    site.location.lon, site.location.lat, site.location.height,
    ra, dec,
    t_mid_0, period, duration,
    max_airmass, minokmoonsep,
    oot_duration, eventclass))

    n_observable = len(event_ind[event_ind])

    # ing_tmid_egr gives time of ingress/tmid/egress. To index with oibeo,
    # repeat it!
    observable_times = ing_tmid_egr[np.repeat(event_ind, 3, axis=0).T]

    lines = []
    for tra_ind in range(n_observable):

        # get ing/mid/egr times, formatted in iso, down to second precision.
        # format them to a width of 22 characters.
        t_ing_str = '{:s}'.format(
            observable_times[tra_ind*3].isot[:-4]).ljust(22,' ').replace('T',' ')
        t_mid_str = '{:s}'.format(
            observable_times[tra_ind*3+1].isot[:-4]).ljust(22,' ').replace('T',' ')
        t_egr_str = '{:s}'.format(
            observable_times[tra_ind*3+2].isot[:-4]).ljust(22,' ').replace('T',' ')

        t_mid_val = observable_times[tra_ind*3+1].jd
        epoch = int(np.round((t_mid_val - t_mid_0) / period.value, 0))

        # if you're not within 1.5 minutes = 1e-3 days of the mid-time being a
        # "true epoch", something is wrong with your calculation. raise error.
        if not (
            np.isclose(
                ((t_mid_val - t_mid_0)/period.value) % 1, 1, atol=1e-3
            ) or
            np.isclose(
                ((t_mid_val - t_mid_0)/period.value) % 1, 0, atol=1e-3
            )
        ):
            raise AssertionError('got error in tra_ind {}')

        epoch_str = '{}'.format(epoch).ljust(10,' ')

        moon_sep_str = '{}d'.format(
            moon_separation[tra_ind]).ljust(5,' ')
        moon_ill_str = '{}%'.format(
            moon_illumination[tra_ind]).ljust(5,' ')

        this_line = '{}{}{}{}{}{}\n'.format(
            epoch_str, t_ing_str, t_mid_str, t_egr_str, moon_sep_str,
            moon_ill_str)

        lines.append(this_line)


    outname = '{}_{}.txt'.format(name, site.name.replace(' ','_'))
    outpath = os.path.join(outdir,outname)
    with open(outpath, mode='w') as f:
        f.writelines(hdr)
        f.writelines(lines)
    print('wrote {}'.format(outpath))

