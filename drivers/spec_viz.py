import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from stringcheese.plotutils import savefig, format_ax

specpath = '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/Veloce/20200130_837.01_Bouma_final_combined.fits'
hdul = fits.open(specpath)

n_orders = hdul[0].shape[0]

for order in range(n_orders):

    start = 200
    end = -200

    flux = hdul[0].data[order, start:end]
    flux_err = hdul[1].data[order, start:end]
    wav = hdul[2].data[order, start:end]

    plt.close('all')
    f,ax = plt.subplots()
    ax.plot(wav, flux, c='k', zorder=3)

    ax.set_xlabel('wavelength [angstrom]')
    ax.set_ylabel('flux [e-]')

    format_ax(ax)
    outpath = '../results/spec_viz/20200130_837.01_order{}.png'.format(str(order).zfill(2))
    savefig(f, outpath)
