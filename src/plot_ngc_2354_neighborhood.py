import os
import numpy as np, matplotlib.pyplot as plt, pandas as pd

from astropy.io.votable import from_table, writeto, parse
from astropy.coordinates import SkyCoord
from astropy import units as u

datadir = '../data/'

def given_votable_get_df(votablepath, assert_equal='source_id'):

    vot = parse(votablepath)
    tab = vot.get_first_table().to_table()
    df = tab.to_pandas()

    if isinstance(assert_equal, str):
        np.testing.assert_array_equal(tab[assert_equal], df[assert_equal])

    return df

nbhrpath = os.path.join(datadir, 'ngc_2354_neighborhood_query-result.vot.gz')
cg18_path = os.path.join(datadir, 'ngc_2354_CG18_subset.vot.gz')

df = given_votable_get_df(nbhrpath, assert_equal='source_id')
cg18_df = given_votable_get_df(cg18_path, assert_equal='Source')

cutoff_probability = 0.1
cg18_df = cg18_df[cg18_df['PMemb'] > cutoff_probability]

target_df = df[df['source_id'] == 5617126180115568256]

##########################################

fig, axs = plt.subplots(2, 3, figsize=(18,12))

######################
# parallax vs pm dec #
######################
ax0 = axs[0,0]

ax0.scatter(
            df['pmdec'],
            df['parallax'],
            c='gray', alpha=0.9, zorder=2, s=15, rasterized=True, linewidths=0,
            label='nbhd stars'
)
ax0.plot(
    target_df['pmdec'],
    target_df['parallax'],
    alpha=1, mew=0.5, zorder=8, label='554-051',
    markerfacecolor='yellow', markersize=20, marker='*', color='black', lw=0
)

ax0.legend(loc='best')
ax0.set_xlabel(r'pmDEC, $\mu_{{\delta}}$ [mas/yr]', fontsize='xx-large')
ax0.set_ylabel('star parallax [mas]', fontsize='xx-large')

ax0.set_xlim([-3,6])

######################
# parallax vs pm ra #
######################
ax1 = axs[0,1]

ax1.scatter(
            df['pmra'],
            df['parallax'],
            c='gray', alpha=0.9, zorder=2, s=15, rasterized=True, linewidths=0,
            label='nbhd stars'
)
ax1.plot(
    target_df['pmra'],
    target_df['parallax'],
    alpha=1, mew=0.5, zorder=8, label='554-051',
    markerfacecolor='yellow', markersize=20, marker='*', color='black', lw=0
)

ax1.legend(loc='best')
ax1.set_xlabel(r'pmRA, $\mu_{{\alpha}} \cos\delta$ [mas/yr]', fontsize='xx-large')
ax1.set_ylabel('star parallax [mas]', fontsize='xx-large')

ax1.set_xlim([-7,3])

##################
# proper motions #
##################
ax2 = axs[0,2]

ax2.scatter(
    df['pmra'],
    df['pmdec'],
    c='gray', alpha=0.9, zorder=2, s=15, rasterized=True, linewidths=0,
    label='nbhd stars'
)
#ax2.scatter(
#    cg18_df['pmRA'],
#    cg18_df['pmDE'],
#    c='k', alpha=0.9, zorder=3, s=15, rasterized=True, linewidths=0,
#    label='CG18 P>{}'.format(cutoff_probability)
#)
ax2.plot(
    target_df['pmra'],
    target_df['pmdec'],
    alpha=1, mew=0.5, zorder=8, label='554-051',
    markerfacecolor='yellow', markersize=20, marker='*', color='black', lw=0
)

ax2.legend(loc='best')

ax2.set_xlabel(r'pmRA, $\mu_{{\alpha}} \cos\delta$ [mas/yr]',
               fontsize='xx-large')
ax2.set_ylabel(r'pmDEC, $\mu_{{\delta}}$ [mas/yr]',
               fontsize='xx-large')

ax2.set_ylim((-20,20))
ax2.set_xlim((-15,15))

ax2.set_ylim((-3,6))
ax2.set_xlim((-7,3))


##############
# HR diagram # 
##############
ax3 = axs[1,0]

ax3.scatter(
    df['phot_bp_mean_mag']-df['phot_rp_mean_mag'], df['phot_g_mean_mag'],
    c='gray', alpha=1., zorder=2, s=15, rasterized=True, linewidths=0,
    label='nbhd stars'
)
ax3.scatter(
    cg18_df['BP-RP'], cg18_df['Gmag'],
    c='k', alpha=1., zorder=3, s=15, rasterized=True, linewidths=0,
    label='CG18 P>{}'.format(cutoff_probability)
)
ax3.plot(
    target_df['phot_bp_mean_mag']-target_df['phot_rp_mean_mag'],
    target_df['phot_g_mean_mag'],
    alpha=1, mew=0.5, zorder=8, label='554-051',
    markerfacecolor='yellow', markersize=20, marker='*', color='black', lw=0
)

ax3.legend(loc='best')
ax3.set_ylabel('G', fontsize='xx-large')
ax3.set_xlabel('Bp - Rp', fontsize='xx-large')

ylim = ax3.get_ylim()
ax3.set_ylim((max(ylim),min(ylim)))
ax3.set_xlim((-0.5, 3.0))

#############
# positions #
#############
ax4 = axs[1,1]

ax4.scatter(
    df['ra'], df['dec'],
    c='gray', alpha=1., zorder=2, s=15, rasterized=True, linewidths=0,
    label='nbhd stars'
)
ax4.scatter(
    cg18_df['RA_ICRS'], cg18_df['DE_ICRS'],
    c='k', alpha=1., zorder=3, s=15, rasterized=True, linewidths=0,
    label='CG18 P>{}'.format(cutoff_probability)
)
ax4.plot(
    target_df['ra'], target_df['dec'],
    alpha=1, mew=0.5, zorder=8, label='554-051',
    markerfacecolor='yellow', markersize=20, marker='*', color='black', lw=0
)

ax4.legend(loc='best')
ax4.set_xlabel(r'RA, $\alpha$ [deg]', fontsize='xx-large')
ax4.set_ylabel('Dec, $\delta$ [deg]', fontsize='xx-large')

# ######################
# # parallax vs DR2 rv #
# ######################
# ax5 = axs[1,2]
# 
# ax5.scatter(
#             df['radial_velocity'],
#             df['parallax'],
#             c='gray', alpha=0.9, zorder=2, s=15, rasterized=True, linewidths=0,
#             label='nbhd stars'
# )
# ax5.plot(
#     30.8,
#     target_df['parallax'],
#     alpha=1, mew=0.5, zorder=8, label='554-051',
#     markerfacecolor='yellow', markersize=20, marker='*', color='black', lw=0
# )
# 
# ax5.legend(loc='best')
# ax5.set_xlabel('DR2 RV (if available) [km/s]', fontsize='xx-large')
# ax5.set_ylabel('star parallax [mas]', fontsize='xx-large')

#######################
# parallax vs RAVE rv #
#######################
ax5 = axs[1,2]

# # query the standard columns and the separation column of RAVE DR5
# from astroquery.vizier import Vizier
# v = Vizier(columns=["*", "+_r"])
# v.ROW_LIMIT = -1
# result = v.query_region(
#         SkyCoord(ra=float(target_df['ra']), dec=float(target_df['dec']),
#                  unit=(u.deg, u.deg), frame='icrs'),
# 	radius=60*u.arcmin,
#         catalog='III/279/rave_dr5'
# )
# res = result[0]
# 
# # have RA, dec in J2000 from RAVE DR5.
# import IPython; IPython.embed()

ax5.scatter(
            df['radial_velocity'],
            df['parallax'],
            c='gray', alpha=0.9, zorder=2, s=15, rasterized=True, linewidths=0,
            label='nbhd stars'
)
ax5.plot(
    30.8,
    target_df['parallax'],
    alpha=1, mew=0.5, zorder=8, label='554-051',
    markerfacecolor='yellow', markersize=20, marker='*', color='black', lw=0
)

ax5.legend(loc='best')
ax5.set_xlabel('DR2 RV (if available) [km/s]', fontsize='xx-large')
ax5.set_ylabel('star parallax [mas]', fontsize='xx-large')



########
# save #
########
for ax in [ax0,ax1,ax2,ax3,ax4,ax5]:
    ax.get_yaxis().set_tick_params(which='both', direction='in',
                                   labelsize='xx-large')
    ax.get_xaxis().set_tick_params(which='both', direction='in',
                                   labelsize='xx-large')

fig.tight_layout()

outpath = '../results/ngc_2354_neighborhood.png'
fig.savefig(outpath, dpi=300, bbox_inches='tight')
print('made {}'.format(outpath))
