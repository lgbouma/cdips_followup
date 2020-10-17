import os, corner, pickle
from datetime import datetime
from glob import glob
import numpy as np, matplotlib.pyplot as plt, pandas as pd, pymc3 as pm
from numpy import array as nparr
from matplotlib import transforms

from aesthetic.plot import savefig, format_ax
from aesthetic.plot import set_style

from cdips.utils.catalogs import get_cdips_catalog
from cdips.utils.gaiaqueries import query_neighborhood

RESULTSDIR = '/Users/luke/Dropbox/proj/cdips_followup/results/TIC_268_neighborhood'

def rainbow_text(x, y, strings, colors, orientation='horizontal',
                 ax=None, **kwargs):
    """
    Take a list of *strings* and *colors* and place them next to each
    other, with text strings[i] being shown in colors[i].

    Parameters
    ----------
    x, y : float
        Text position in data coordinates.
    strings : list of str
        The strings to draw.
    colors : list of color
        The colors to use.
    orientation : {'horizontal', 'vertical'}
    ax : Axes, optional
        The Axes to draw into. If None, the current axes will be used.
    **kwargs
        All other keyword arguments are passed to plt.text(), so you can
        set the font size, family, etc.
    """
    if ax is None:
        ax = plt.gca()
    t = ax.transAxes
    canvas = ax.figure.canvas

    props = dict(boxstyle='square', facecolor='white', alpha=1, pad=0.05,
                 linewidth=0)

    kwargs.update(rotation=0, va='bottom', ha='right')

    for s, c in zip(strings, colors):
        text = ax.text(x, y, s + " ", color=c, transform=t, bbox=props,
                       **kwargs)

        # Need to draw to update the text position.
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()
        t = transforms.offset_copy(
            text.get_transform(), y=0.7*ex.height, units='dots'
        )


def plot_TIC268_nbhd_small(outdir=RESULTSDIR):

    set_style()

    df = get_cdips_catalog(ver=0.4)

    kc19_sel = (
        (df.cluster.str.contains('NGC_2516')) &
        (df.reference.str.contains('Kounkel_2019'))
    )
    cg18_sel = (
        (df.cluster.str.contains('NGC_2516')) &
        (df.reference.str.contains('CantatGaudin_2018'))
    )

    target_df = df[df.source_id == 5489726768531119616] # TIC 2683...
    kc19_df = df[kc19_sel]
    cg18_df = df[cg18_sel]

    kc19_df = kc19_df[~(kc19_df.source_id.isin(cg18_df.source_id))]

    ##########

    # NGC 2516 rough
    bounds = {
        'parallax_lower': 1.5,
        'parallax_upper': 4.0,
        'ra_lower': 108,
        'ra_upper': 132,
        'dec_lower': -76,
        'dec_upper': -45
    }
    groupname = 'customngc2516'

    nbhd_df = query_neighborhood(bounds, groupname, n_max=6000,
                                 overwrite=False, manual_gmag_limit=17)

    sel_nbhd = (
        (~nbhd_df.source_id.isin(kc19_df.source_id))
        &
        (~nbhd_df.source_id.isin(cg18_df.source_id))
    )
    from copy import deepcopy
    orig_nbhd_df = deepcopy(nbhd_df)
    nbhd_df = nbhd_df[sel_nbhd]

    print(f'Got {len(nbhd_df)} neighbors')
    print(f'Got {len(cg18_df)} in core')
    print(f'Got {len(kc19_df)} in corona')

    ##########

    plt.close('all')

    f, axs = plt.subplots(figsize=(4,3), ncols=2)

    xv, yv = 'ra', 'dec'
    axs[0].scatter(
        nbhd_df[xv], nbhd_df[yv], c='gray', alpha=0.9, zorder=2, s=7,
        rasterized=True, linewidths=0, label='Field', marker='.'
    )
    axs[0].scatter(
        kc19_df[xv], kc19_df[yv], c='C0', alpha=0.9, zorder=3, s=7,
        rasterized=True, linewidths=0, label='Corona', marker='.'
    )
    axs[0].scatter(
        cg18_df[xv], cg18_df[yv], c='k', alpha=0.9, zorder=4, s=7,
        rasterized=True, label='Core', marker='.'
    )
    axs[0].plot(
        target_df[xv], target_df[yv], alpha=1, mew=0.5,
        zorder=8, label='TOI 1937', markerfacecolor='lightskyblue',
        markersize=14, marker='*', color='black', lw=0
    )

    axs[0].set_xlabel(r'$\alpha$ [deg]')
    axs[0].set_ylabel(r'$\delta$ [deg]')
    axs[0].set_xlim([108, 132])
    axs[0].set_ylim([-76, -45])

    ##########

    get_yval = (
        lambda _df: np.array(
            _df['phot_g_mean_mag'] + 5*np.log10(_df['parallax']/1e3) + 5
        )
    )
    get_xval = (
        lambda _df: np.array(
            _df['phot_bp_mean_mag'] - _df['phot_rp_mean_mag']
        )
    )

    axs[1].scatter(
        get_xval(nbhd_df), get_yval(nbhd_df), c='gray', alpha=0.9, zorder=2,
        s=7, rasterized=True, linewidths=0, label='Field', marker='.'
    )
    axs[1].scatter(
        get_xval(kc19_df), get_yval(kc19_df), c='C0', alpha=1, zorder=3,
        s=7, rasterized=True, linewidths=0, label='Corona', marker='.'
    )
    axs[1].scatter(
        get_xval(cg18_df), get_yval(cg18_df), c='k', alpha=0.9,
        zorder=4, s=7, rasterized=True, linewidths=0, label='Core', marker='.'
    )
    axs[1].plot(
        get_xval(target_df), get_yval(target_df), alpha=1, mew=0.5,
        zorder=8, label='TOI 1937', markerfacecolor='lightskyblue',
        markersize=14, marker='*', color='black', lw=0
    )

    axs[1].set_ylim(axs[1].get_ylim()[::-1])

    axs[1].set_xlabel('Bp - Rp [mag]')
    axs[1].set_ylabel('Absolute G [mag]', labelpad=-6)

    ##########

    words = ['Field', 'Corona', 'Core', 'TOI1937'][::-1]
    colors = ['gray', 'C0', 'k', 'lightskyblue'][::-1]
    rainbow_text(0.98, 0.02, words, colors, size='medium', ax=axs[0])

    f.tight_layout(w_pad=2)

    outpath = os.path.join(outdir, 'small_ngc2516.png')
    savefig(f, outpath)


if __name__=="__main__":
    plot_TIC268_nbhd_small()
    #plot_TIC268_nbhd_full()
