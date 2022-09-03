import numpy as np
from scipy.optimize import brentq

iodfactor = 0.7

def exposure_time(vmag, counts, iod=False, time1=110.0, vmag1=8.0, exp1=250.0):
    """Exposure Time

    Estimate exposure time based on scaling. Cannonical exposure time
    is 110s to get to 250k on 8th mag star with iodine cell in.

    Args:
        vmag (float): V-band magnitude
            250 = 250k, 10 = 10k (CKS) i.e. SNR = 45 per pixel.
        counts (float): desired number of counts.
            250 = 250k, 10 = 10k (CKS) i.e. SNR = 45 per pixel.
        iod (bool): is iodine cell in or out? If out, throughput is higher
            by 30%

    Returns:
        float: exposure time (seconds)

    """

    # flux star / flux 8th mag star
    fluxfactor = 10.0**(-0.4*(vmag-vmag1))
    time = time1 / fluxfactor
    time *= counts / exp1
    if iod==False:
        time *= iodfactor
    return time

def exposure_counts(vmag, time, **kwargs):
    """Exposure counts

    Inverse of `exposure_time.` Given a magnitude and an exposure
    time, how many counts will be collected?

    Args:
        vmag (float) : vband magnitude
        time (float) : exposure time (seconds)
        **kwargs : keyword arguments passed to exposure_time

    Returns:
        float: expected number of counts

    """
    f = lambda counts : exposure_time(vmag, counts, **kwargs) - time
    _counts = brentq(f,0,2000,)
    return _counts
