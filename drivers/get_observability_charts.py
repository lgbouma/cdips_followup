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


def get_observability_fraction(name="WASP 4", site='keck', ra=None, dec=None,
                               start_time=Time('2019-09-13 20:00:00'),
                               end_time=Time('2020-07-31 20:00:00')):

    if isinstance(name,str) and ra is None and dec is None:
        target = FixedTarget.from_name(name)
    elif isinstance(ra,float) and isinstance(dec,float):
        target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        target = FixedTarget(coord=target_coord, name=name)
    else:
        raise NotImplementedError('failed to make target')

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

    return frac, ever_observable, days



def make_observability_chart_singlesite(
    name="WASP 4", site='keck', ra=None, dec=None,
    start_time=Time('2019-09-13 20:00:00'),
    end_time=Time('2020-07-31 20:00:00'),
    outdir=None):

    frac, ever_observable, days = get_observability_fraction(
        name=name, site=site, ra=ra, dec=dec,
        start_time=start_time, end_time=end_time
    )

    #
    # make the plot
    #
    plt.close('all')
    fig, ax = plt.subplots(figsize=(4,3))

    ax.plot_date(days.plot_date, frac*24, label=name, lw=1,
                 alpha=1, fmt='k-')

    ax.set_xlim([days[0].plot_date, days[-1].plot_date])
    ax.set_title(
        '{}: altitude 20-85deg, airmass<3, w/in civil twilight'.format(name),
        fontsize='xx-small'
    )
    plt.setp(ax.get_xticklabels(), rotation=90, ha='right')

    ax.set_xlabel('time')
    ax.set_ylabel('hours visible per night')
    ax.set_ylim((0,12))

    outpath = os.path.join(
        outdir,
        '{:s}_frac_visible_{:s}.png'.format(name, site.replace(' ','_').lower())
    )
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print('saved {}'.format(outpath))


def make_observability_chart_multisite(
    name="WASP 4", sites=['keck','las campanas'], ra=None, dec=None,
    start_time=Time('2019-09-13 20:00:00'),
    end_time=Time('2020-07-31 20:00:00'),
    outdir=None,
    overwrite=False,
    save_csv=True):

    outpath = os.path.join(
        outdir,
        '{}_frac_visible_multisite.png'.format(name)
    )

    if os.path.exists(outpath) and not overwrite:
        print('found {}, skipping'.format(outpath))
        return

    fracs = []
    for site in sites:
        frac, _, days = get_observability_fraction(
            name=name, site=site, ra=ra, dec=dec,
            start_time=start_time, end_time=end_time
        )
        fracs.append(frac)

    radecstr = 'RA {:.1f}, dec {:.1f}'.format(float(ra), float(dec))

    #
    # make the plot
    #
    plt.close('all')
    fig, ax = plt.subplots(figsize=(4,3))

    for ix, frac in enumerate(fracs):
        ax.plot_date(days.plot_date, frac*24, lw=1, alpha=1,
                     fmt='C{}-'.format(ix), label=sites[ix])

    ax.set_xlim([days[0].plot_date, days[-1].plot_date])
    ax.set_title(
        '{}: alt 20-85deg, airmass<3, civil twilight'.format(name),
        fontsize='xx-small'
    )
    ax.text(0.03, 0.97, radecstr, va='top', ha='left', transform=ax.transAxes)
    plt.setp(ax.get_xticklabels(), rotation=75, ha='right', fontsize='small')
    plt.setp(ax.get_yticklabels(), fontsize='small')

    ax.legend(loc='upper right', fontsize='xx-small')

    ax.set_xlabel('time')
    ax.set_ylabel('hours visible per night')
    ax.set_ylim((0,12))

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print('saved {}'.format(outpath))

    if save_csv:
        outcsv = outpath.replace('.png','.csv')
        outdf = pd.DataFrame({
            'plot_date':days.plot_date,
        })
        for ix, frac in enumerate(fracs):
            label = sites[ix]
            outdf[label] = frac*24
        outdf.to_csv(outcsv, index=False)
        print('made {}'.format(outcsv))


if __name__ == "__main__":

    name='WASP 4'

    make_airmass_chart(
        name=name, site='keck',
        ra=None, dec=None,
        outdir='../results/followup_planning/WASP_4',
        start_time=Time('2019-09-13 20:00:00', format='iso'),
        end_time=Time('2020-09-13 20:00:00', format='iso')
    )
