import os, corner, pickle
from datetime import datetime
from glob import glob
import numpy as np, matplotlib.pyplot as plt, pandas as pd, pymc3 as pm
from numpy import array as nparr

from aesthetic.plot import savefig, format_ax
from aesthetic.plot import set_style

from cdips.utils.catalogs import get_cdips_catalog

RESULTSDIR = '/Users/luke/Dropbox/proj/cdips_followup/results/TIC_268_neighborhood'

def plot_TIC268_nbhd_full(outdir=RESULTSDIR):

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

    ##########

    plt.close('all')

    params = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']
    nparams = len(params)

    qlimd = {
        'ra': 0,
        'dec': 0,
        'parallax': 0,
        'pmra': 0,
        'pmdec': 0,
        'radial_velocity': 0
    } # whether to limit axis by quantile

    ldict = {
        'ra': r'$\alpha$ [deg]',
        'dec': r'$\delta$ [deg]',
        'parallax': r'$\pi$ [mas]',
        'pmra': r'$\mu_{{\alpha}} \cos\delta$ [mas/yr]',
        'pmdec':  r'$\mu_{{\delta}}$ [mas/yr]',
        'radial_velocity': 'RV [km/s]'
    }

    f, axs = plt.subplots(figsize=(6,6), nrows=nparams-1, ncols=nparams-1)

    for i in range(nparams):
        for j in range(nparams):
            print(i,j)
            if j == nparams-1 or i == nparams-1:
                continue
            if j>i:
                axs[i,j].set_axis_off()
                continue

            xv = params[j]
            yv = params[i+1]
            print(i,j,xv,yv)

            axs[i,j].scatter(
                kc19_df[xv], kc19_df[yv], c='gray', alpha=0.9, zorder=2, s=5,
                rasterized=True, linewidths=0, label='Core', marker='.'
            )
            axs[i,j].scatter(
                cg18_df[xv], cg18_df[yv], c='k', alpha=0.9,
                zorder=3, s=5, rasterized=True, linewidths=0, label='Corona'
            )
            axs[i,j].plot(
                target_df[xv], target_df[yv], alpha=1, mew=0.5,
                zorder=8, label='TOI 1937', markerfacecolor='yellow',
                markersize=14, marker='*', color='black', lw=0
            )

            # set the axis limits as needed
            if qlimd[xv]:
                xlim = (np.nanpercentile(kc19_df[xv], 25),
                        np.nanpercentile(kc19_df[xv], 75))
                axs[i,j].set_xlim(xlim)
            if qlimd[yv]:
                ylim = (np.nanpercentile(kc19_df[yv], 25),
                        np.nanpercentile(kc19_df[yv], 75))
                axs[i,j].set_ylim(ylim)

            # fix labels
            if j == 0 :
                axs[i,j].set_ylabel(ldict[yv])
                if not i == nparams - 2:
                    # hide xtick labels
                    labels = [item.get_text() for item in axs[i,j].get_xticklabels()]
                    empty_string_labels = ['']*len(labels)
                    axs[i,j].set_xticklabels(empty_string_labels)

            if i == nparams - 2:
                axs[i,j].set_xlabel(ldict[xv])
                if not j == 0:
                    # hide ytick labels
                    labels = [item.get_text() for item in axs[i,j].get_yticklabels()]
                    empty_string_labels = ['']*len(labels)
                    axs[i,j].set_yticklabels(empty_string_labels)

            if (not (j == 0)) and (not (i == nparams - 2)):
                # hide ytick labels
                labels = [item.get_text() for item in axs[i,j].get_yticklabels()]
                empty_string_labels = ['']*len(labels)
                axs[i,j].set_yticklabels(empty_string_labels)
                # hide xtick labels
                labels = [item.get_text() for item in axs[i,j].get_xticklabels()]
                empty_string_labels = ['']*len(labels)
                axs[i,j].set_xticklabels(empty_string_labels)

    # axs[2,2].legend(loc='best', handletextpad=0.1, fontsize='medium', framealpha=0.7)
    # leg = axs[2,2].legend(bbox_to_anchor=(0.8,0.8), loc="upper right",
    #                       handletextpad=0.1, fontsize='medium',
    #                       bbox_transform=f.transFigure)


    for ax in axs.flatten():
        format_ax(ax)

    f.tight_layout(h_pad=0.1, w_pad=0.1)

    outpath = os.path.join(outdir, 'full_kinematics.png')
    savefig(f, outpath)



