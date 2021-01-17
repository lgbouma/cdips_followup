import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
from glob import glob

from aesthetic.plot import savefig, format_ax, set_style

def plot_dilution_histogram(ap=None):
    '''
    Make the plot of log10(dilution [2px aperture]) vs log10(distance [pc]) for
    T<16 mag, d<2kpc cluster members.
    '''

    df = pd.read_csv('../../data/cg18_cdips_table1_subset_with_dilution.csv')

    sel = (df.parallax > 0.5) & (df.Tmag_pred < 16) # d < 2kpc, T<16
    df = df[sel]

    dil_n_px = np.array(df['dilution_ap{:.2f}'.format(ap)]) # y
    dil_n_px[dil_n_px > 0.999 ] = 0.999

    set_style()
    plt.close('all')
    fig, ax = plt.subplots(figsize=(3,3))

    ax.set_xscale('log')
    ax.set_xlabel(
        '(Target flux)/(Total flux in aperture)'.
        format(ap)
    )
    ax.set_ylabel('Relative probability')

    ax.set_xlim((10**(-2.05), 1.1))
    ax.tick_params(which='both', direction='in', zorder=0)

    xmin, xmax = 10**(-3), 10**1

    log_dil_n_px_bins = np.linspace(np.log10(xmin), np.log10(xmax), 17)

    x = 10**log_dil_n_px_bins
    y = np.histogram(np.log10(dil_n_px), log_dil_n_px_bins)[0]
    x = np.array(list(zip(x[:-1], x[1:]))).flatten()
    y = np.array(list(zip(y, y))).flatten()
    ax.plot(x, y, lw=1, color='black')
    inds = (x <= 0.1)
    ax.fill_between(x[inds], y[inds], np.zeros_like(y[inds]), facecolor='none',
            hatch='/', edgecolor='gray', lw=0)

    frac_not_ok = np.sum(y[inds]) / np.sum(y)
    nonrecov_str = r'$\approx$'+'{:d}%\ntoo crowded'.format(int(100*frac_not_ok))
    recov_str = r'$\approx$'+'{:d}%\nrecoverable'.format(
            int(round(100*(1-frac_not_ok))))

    t = ax.text(10**(-0.5), 11500, recov_str,
            verticalalignment='center',horizontalalignment='center',fontsize='medium')
    t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='gray'))
    t= ax.text(10**(-1.5), 11500, nonrecov_str,
            verticalalignment='center',horizontalalignment='center',fontsize='medium')
    t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='gray'))

    #ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))

    format_ax(ax)
    # ax.yaxis.set_ticks_position('both')
    # ax.xaxis.set_ticks_position('both')
    # ax.get_yaxis().set_tick_params(which='both', direction='in')
    # ax.get_xaxis().set_tick_params(which='both', direction='in')
    # ax.xaxis.set_tick_params(labelsize='small')
    # ax.yaxis.set_tick_params(labelsize='small')

    # ax.tick_params(which='both', direction='in', zorder=0)
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)

    ax.set_ylim((0, max(ax.get_ylim())))

    outname = (
        '../../results/dilution_histogram/'
        'dil_Tmaglt16_dlt2kpc_cg18_ap{:.2f}px.png'.
        format(ap)
    )

    savefig(fig, outname)

if __name__ == "__main__":

    aps = np.arange(0.75, 2.75, 0.25)
    for ap in aps:
        plot_dilution_histogram(ap=ap)
