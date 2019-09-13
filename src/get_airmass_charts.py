import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
from glob import glob

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

from astroplan.plots import plot_airmass
from astroplan import (FixedTarget, Observer, months_observable,
                       AtNightConstraint, AltitudeConstraint)

import datetime as dt

def make_airmass_chart(name="WASP 4", site='keck', ra=None, dec=None,
                       start_time=Time('2019-09-13 20:00:00'),
                       end_time=Time('2020-07-31 20:00:00'),
                       outdir='../results/followup_planning/WASP_4'):
    """
    2020A:
        - Feb 1 - Jul 31 2020
    """

    if (
        not (isinstance(name,str)
             and ra is None
             and dec is None)
        and not (isinstance(name,None)
             and isinstance(ra,float)
             and isinstance(dec,float))
    ):
        raise NotImplementedError('Either give name, or give ra and dec.')

    if isinstance(name, str):
        target = FixedTarget.from_name(name)
    elif isinstance(ra,float) and isinstance(dec,float):
        #FIXME
        target = FixedTarget.from_name(name)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    observer = Observer.at_site(site)

    constraints = [AltitudeConstraint(min=20*u.deg, max=85*u.deg),
                   AtNightConstraint.twilight_civil()]

    observe_time = Time(
        np.arange(start_time.decimalyear, end_time.decimalyear,
                  0.5/(24*365.25)),
        format='decimalyear'
    )

    #FIXME this isn't working ?!
    import IPython; IPython.embed()

    # FIXME want this?
    # computed observability on "best_months" grid of 0.5 hr
    # best_months = months_observable(constraints, observer, [target])
    # print('for {}, got best-months on 0.5 hour grid:'.format(name))
    # print(best_months)
    # print('where 1 = Jan, 2 = Feb, etc.')

    plt.close('all')

    plot_airmass(target, observer, observe_time)

    outpath = os.path.join(outdir, 'airmass.png')
    plt.savefig(outpath, dpi=300)
    print('saved {}'.format(outpath))



if __name__ == "__main__":

    name='WASP 4'

    make_airmass_chart(name=name, outdir='../results/followup_planning/WASP_4',
                       start_time=Time('2018-09-13 20:00:00', format='iso'),
                       end_time=Time('2021-09-13 20:00:00', format='iso')
                      )
