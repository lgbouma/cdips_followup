import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
from glob import glob

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

from astroplan import (
    FixedTarget, Observer, months_observable, AtNightConstraint,
    AltitudeConstraint, is_observable, observability_table,
    AirmassConstraint
)

def make_airmass_chart(name="WASP 4", site='keck', ra=None, dec=None,
                       start_time=Time('2019-09-13 20:00:00'),
                       end_time=Time('2020-07-31 20:00:00'),
                       outdir='../results/followup_planning/WASP_4',
                       check_months_observable=True):
    """
    2020A:
        - Feb 1 - Jul 31 2020
    """

    if isinstance(name,str) and ra is None and dec is None:
        target = FixedTarget.from_name(name)
    elif isinstance(ra,float) and isinstance(dec,float):
        target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        target = FixedTarget(coord=target_coord, name=name)
    else:
        raise NotImplementedError('failed to make target')

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    observer = Observer.at_site(site)

    constraints = [AltitudeConstraint(min=20*u.deg, max=85*u.deg),
                   AirmassConstraint(3),
                   AtNightConstraint.twilight_civil()]

    # over every day between start and end time, check if the observing
    # constraints are meetable.
    days = Time(
        np.arange(start_time.decimalyear, end_time.decimalyear,
                  1/(365.25)),
        format='decimalyear'
    )

    frac, ever_observable = [], []

    for day in days:

        table = observability_table(constraints, observer, [target],
                                    time_range=day)
        frac.append(float(table['fraction of time observable']))
        ever_observable.append(bool(table['ever observable']))

    ever_observable = np.array(ever_observable)
    frac = np.array(frac)

    #
    # make the plot
    #
    plt.close('all')
    fig, ax = plt.subplots(figsize=(4,3))

    ax.plot_date(days.plot_date, frac, label=name, lw=1,
                 alpha=1, fmt='k-')

    ax.set_xlim([days[0].plot_date, days[-1].plot_date])
    ax.set_title(
        '{}: altitude 20-85deg, airmass<3, w/in civil twilight'.format(name),
        fontsize='small'
    )
    plt.setp(ax.get_xticklabels(), rotation=90, ha='right')

    ax.set_xlabel('time')
    ax.set_ylabel('fraction visible')

    outpath = os.path.join(
        outdir,
        'frac_visible_{:s}.png'.format(site.replace(' ','_').lower())
    )
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print('saved {}'.format(outpath))

    #
    # optionally computes observability on "best_months" grid of 0.5 hr
    #
    if check_months_observable:
        best_months = months_observable(constraints, observer, [target])
        outfile = os.path.join(outdir,'months_observable.txt')
        with open(outfile, 'w') as f:
            f.write('for {}, got best-months on 0.5 hour grid:'.format(name))
            f.write(repr(best_months))
            f.write('where 1 = Jan, 2 = Feb, 3=March, etc.')
        print('made {}'.format(outfile))


if __name__ == "__main__":

    name='WASP 4'

    make_airmass_chart(
        name=name, site='keck',
        ra=None, dec=None,
        outdir='../results/followup_planning/WASP_4',
        start_time=Time('2019-09-13 20:00:00', format='iso'),
        end_time=Time('2020-09-13 20:00:00', format='iso'),
        check_months_observable=False
    )
