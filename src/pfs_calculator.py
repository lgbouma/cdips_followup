#
# for an 8m/s internal precision measurement, according to Johanna's empirical
# plots
#
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from numpy import array as nparr

def log_interp1d(x, y, kind='linear'):
    logx = np.log10(x)
    logy = np.log10(y)
    lin_interp = interp1d(logx, logy, kind=kind)
    log_interp = lambda z: np.power(10., lin_interp(np.log10(z)))
    return log_interp

def given_Vmag_get_PFS_exptime_minutes(Vmag):

    df = pd.DataFrame({
        'Vmag': [8, 10, 12, 12.5, 15],
        'exptime_min': [0.15, 0.8, 5, 8, 80]
    })

    fn = log_interp1d(nparr(df['Vmag']), nparr(df['exptime_min']))

    return fn(Vmag)


def test_PFS_calculator():

    _Vmag = np.linspace(8,15,100)

    exptime = given_Vmag_get_PFS_exptime_minutes(_Vmag)

    f,ax = plt.subplots(figsize=(4,3))
    ax.plot(_Vmag, exptime)

    ax.set_xlabel('Vmag')
    ax.set_ylabel('exptime [min]')
    ax.set_yscale('log')

    f.savefig('../results/test_PFS_calculator.png', bbox_inches='tight', dpi=250)

if __name__=="__main__":
    test_PFS_calculator()
