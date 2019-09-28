# -*- coding: utf-8 -*-

import numpy as np

def get_chiron_reqd_exptime(gmag,
                            mode,
                            target_snr_per_resoln_element_all_exposures=50,
                            n_exposures=1,
                            verbose=True):
    """
    calculate time per exposure needed to get one high-resolution spectrum
    (i.e.  one "good" RV point), when binning together n_exposures.
    """
    # assume: binning all exposures to get one high-resolution spectrum.

    assert mode == 'slicer_mode' or mode == 'fiber_mode'

    slicer_resolution = 8e4
    fiber_resolution = 2.8e4

    if mode == 'slicer_mode':
        # N pixels per resolution element = N_pixels / (resolving power * N_orders).
	# Where N_pixels = 1.6e7
	# resolving power = R = λ/Δλ = 8e4
	# N_orders ~= 70.
	# 
	# Then
	# 
	#  N_pixels / (resolving power * N_orders) ~= 1.6e7/(8e4*70) ~= 3.
        pixels_per_res_element = 3

    elif mode == 'fiber_mode':
        # if you go to lower resolution, of 28k instead of 80k, you get
        # linearly more pixels per resolution element.
        pixels_per_res_element = 3*slicer_resolution/fiber_resolution

    target_snr_per_pixel_per_exposure = (
        target_snr_per_resoln_element_all_exposures /
        ( np.sqrt(pixels_per_res_element * n_exposures) )
    )
    if verbose:
        print('G={}, mode: {}'.format(gmag, mode))
        print('target snr per resolution element per exposure: {:.1f}'.format(
            target_snr_per_resoln_element_all_exposures/np.sqrt(n_exposures)))
        print('target snr per pixel per exposure: {:s}'.format(
            repr(target_snr_per_pixel_per_exposure)))

    # photons/s/pixel is flux of V=0 star per 0.0202Angstrom pixel. the number
    # of angstroms per pixel is invariant to the spectograph resolution.
    F_0 = 3.4e5
    # note: 5500Angstrom / 8e4 resolution = 0.06875 angstroms, i.e. F_0 is the
    # flux of a V=0 per LITERAL PIXEL (not resolution element) in slicer mode.
    # There are ~3 pixels per resolution element in slicer mode.
    #
    # In the fiber mode, 5500Angstrom / 2.8e4 resolution = 0196Angstroms -- the
    # same number of photons per second are hitting each literal pixel. However
    # the ability to distinguish them in resolution elements has worsened.

    if mode == 'slicer_mode':
        R = 8 # CCD readout noise in electrons
        K = 9 # number of binned pixels across the order
        efficiency = 0.06 # total efficiency
    elif mode == 'fiber_mode':
        R = 4.3 # CCD readout noise in electrons
        K = 2.5 # number of binned pixels across the order
        efficiency = 0.06*(1/.82) # total efficiency, Table2 of Tokovinin+2013

        # FIXME
        # i don't understand why, but everything in fiber mode is a ~factor of 2
        # off from the plots in the manual otherwise :-(
        fudge = 2
        efficiency *= fudge
        # FIXME

    # Gamma = flux = count rate = photons/sec/px
    Gamma = F_0 * efficiency * 10**(-0.4*gmag)

    snr = target_snr_per_pixel_per_exposure

    # Take Equation 1 of Tokovinin et al, note that N_ph = Gamma *
    # time_per_exposure, then solve the quadratic equation for
    # time_per_exposure.
    # You get the following.
    time_per_exposure = 0.5*(
        snr**2/Gamma
        +
        np.sqrt( (snr**2/Gamma)**2 + 4*snr**2*K*R**2/(Gamma**2) )
    )

    return time_per_exposure
