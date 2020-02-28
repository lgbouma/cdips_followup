from math import pi
import numpy as np, matplotlib.pyplot as plt, pandas as pd
import os

from astropy.time import Time

from astropy.time import Time
from astropy.coordinates import get_body, get_sun, get_moon, SkyCoord
import astropy.units as u

from astroplan import (FixedTarget, Observer, EclipsingSystem,
                       PrimaryEclipseConstraint, is_event_observable,
                       AtNightConstraint, AltitudeConstraint,
                       LocalTimeConstraint, MoonSeparationConstraint,
                       AirmassConstraint, moon)
from astroplan.plots import plot_airmass

from matplotlib.dates import (YEARLY, DateFormatter, rrulewrapper,
                              RRuleLocator, drange)
from matplotlib.ticker import AutoMinorLocator

import datetime as dt



def main():

    ticid = '268301217'
    ra = '07:45:28.97'
    dec = '-52:22:59.73'

    t0 = 2458860.60021*u.day # BJD_TDB
    period = 0.94667643*u.day

    site = Observer.at_site('Cerro Tololo')
    start_time = Time('2020-03-02 22:00:00')
    end_time = Time('2020-03-17 07:00:00')

    # axvspan windows
    obs_per_night = [(
        Time('2020-03-{}'.format(str(d).zfill(2)) + ' 00:00:00'),
        Time('2020-03-{}'.format(str(d).zfill(2)) + ' 04:00:00')
        )
        for d in range(2, 18)
    ]

    plot_quadrature_windows(ra, dec, ticid, t0, period, start_time, end_time,
                            site, obs_per_night=obs_per_night, N_points=500)


def plot_quadrature_windows(ra, dec, ticid, t0, period, start_time, end_time, site,
                            obs_per_night=None, N_points=500):

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

    target = FixedTarget(coord=target_coord, name='TIC'+ticid)

    delta_t = end_time - start_time
    observe_time = start_time + delta_t*np.linspace(0, 1, N_points)

    fig = plt.figure(figsize=(24,4))
    ax = plt.gca()

    plot_airmass(target, site, observe_time)

    plt.plot_date(observe_time.plot_date, np.ones(len(observe_time)), alpha=0)

    formatter = DateFormatter('%m/%d/%y')
    ax.xaxis.set_major_formatter(formatter)

    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylabel('Airmass', color='C0')

    if not obs_per_night is None:

        ix = 0
        for o in obs_per_night:

            if ix == 0:
                ax.plot_date([o[0].plot_date, o[1].plot_date], [1.01, 1.01],
                             c='C2', alpha=0.5, zorder=2, fmt='-', lw=10,
                             label='00 to 04 UT')
            else:
                ax.plot_date([o[0].plot_date, o[1].plot_date], [1.01, 1.01],
                             c='C2', alpha=0.5, zorder=2, fmt='-', lw=10)

            ix += 1

    # calculate the expected RV signal. note we are ignoring barycentric timing
    # effects (BJD vs JD), which are of order 16 minutes.
    rv = np.sin( ( (2*pi/period) * (observe_time.jd*u.day - t0) )*u.rad )

    ax.legend(loc='best')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax2 = ax.twinx()
    ax2.plot_date(observe_time.plot_date, rv, fmt='-', c='C1')
    ax2.set_ylabel('RV', color='C1')

    outpath = (
        '../../results/followup_planning/202002_pfs_quadrature/TIC{}.png'.
        format(ticid)
    )
    fig.savefig(outpath, bbox_inches='tight')
    print('made {}'.format(outpath))


if __name__ == "__main__":
    main()
