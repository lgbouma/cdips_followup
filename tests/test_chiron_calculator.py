from cdips_followup.exposure_calculators.chiron_calculator import get_chiron_reqd_exptime

Gmag = 7.5
mode='slicer_mode'
exptime = get_chiron_reqd_exptime(Gmag,
                                  mode,
                                  target_snr_per_resoln_element_all_exposures=50,
                                  n_exposures=1,
                                  verbose=True)

# in 1 minute, slicer mode gives SNR of 50 across the resolution elements
# SNR of 50/sqrt(3) ~= 28.9, which is over ~3 pixels (per resolution element).

print(exptime)

# this matches Fig 10 of Tokovinin+2013
assert exptime - 60 < 0.1

#==========================================
Gmag = 10.0
mode='slicer_mode'
exptime = get_chiron_reqd_exptime(Gmag,
                                  mode,
                                  target_snr_per_resoln_element_all_exposures=50,
                                  n_exposures=1,
                                  verbose=True)
print(exptime)
#==========================================
Gmag = 10.0
mode='fiber_mode'
exptime = get_chiron_reqd_exptime(Gmag,
                                  mode,
                                  target_snr_per_resoln_element_all_exposures=50,
                                  n_exposures=1,
                                  verbose=True)
print(exptime)
#==========================================
Gmag = 12.5
mode='fiber_mode'
exptime = get_chiron_reqd_exptime(Gmag,
                                  mode,
                                  target_snr_per_resoln_element_all_exposures=50,
                                  n_exposures=1,
                                  verbose=True)
print(exptime)

#==========================================
Gmag = 14.0
mode='fiber_mode'
exptime = get_chiron_reqd_exptime(Gmag,
                                  mode,
                                  target_snr_per_resoln_element_all_exposures=50,
                                  n_exposures=1,
                                  verbose=True)
print(exptime)

#==========================================
Gmag = 15.0
mode='fiber_mode'
exptime = get_chiron_reqd_exptime(Gmag,
                                  mode,
                                  target_snr_per_resoln_element_all_exposures=50,
                                  n_exposures=1,
                                  verbose=True)
print(exptime)
