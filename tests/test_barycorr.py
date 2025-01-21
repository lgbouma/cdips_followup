from __future__ import print_function
from __future__ import division
from astropy.time import Time
import numpy as np
from barycorrpy import get_BC_vel , exposure_meter_BC_vel
#from . import utc_tdb

def run_sample():
        a=[]
        b=0

        JDUTC = 2458000 # Also accepts float input for JDUTC. Verify scale and format

        # Observation of Tau Ceti taken from CTIO on JD 2458000.
        # Observatory location manually entered. Stellar positional parameters taken from Hipparcos Catalogue
        result = get_BC_vel(JDUTC=JDUTC, hip_id=8102, lat=-30.169283,
                            longi=-70.806789, alt=2241.9, ephemeris='de430',
                            zmeas=0.0)

        if np.isclose(a = result[0], b = 15403.9508, atol = 1e-2, rtol = 0):
            a.append('result')
            b+=1


        # Observation of Tau Ceti taken from CTIO on JD 2458000. Observatory location taken from Astropy list.
        # Stellar positional parameters taken from Hipparcos Catalogue
        JDUTC = Time(2458000, format='jd', scale='utc')
        result2  = get_BC_vel(JDUTC=JDUTC, hip_id=8102, obsname='CTIO',
                              ephemeris='de430')

        if np.isclose(a = result2[0], b = 15403.9608, atol = 1e-2, rtol = 0):
            a.append('result2')
            b+=1

        # Observation of Tau Ceti taken from CTIO on JD 2458000,2458010,2458020.
        # Observatory and stellar parameters entered by user.
        # Use DE405 ephemeris

        obsname=''
        lat=-30.169283
        longi=-70.806789
        alt=2241.9

        epoch = 2451545.0

        ra=26.0213645867
        dec=-15.9395557246
        pmra = -1721.05
        pmdec = 854.16
        px = 273.96

        rv = 0.0
        zmeas=0.0
        # Can also enter JDUTC as float instead of Astropy Time Object
        JDUTC=[2458000,2458000.00001,2458000.00002]
        ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp'

        result3=get_BC_vel(JDUTC=JDUTC, ra=ra, dec=dec, obsname=obsname,
                           lat=lat, longi=longi, alt=alt, pmra=pmra,
                           pmdec=pmdec, px=px, rv=rv, zmeas=zmeas,epoch=epoch,
                           ephemeris=ephemeris, leap_update=True)

        if np.allclose([result3[0][0],result3[0][1],result3[0][2]],
                       [15407.4860,15407.4723,15407.4586],atol = 1e-2, rtol = 0):

            a.append('result3')
            b+=1

        # Predictive Mode
        result5 = get_BC_vel(JDUTC=2458000, hip_id=8102, lat=-30.169283,
                             longi=-70.806789, alt=2241.9, ephemeris='de430',
                             zmeas=0.0, predictive=True)

        if np.isclose(a = result5[0], b = -15403.15938, atol = 1e-2, rtol = 0):
            a.append('result5')
            b+=1


        import IPython; IPython.embed()

        result6 = get_BC_vel(JDUTC=2458000, lat=-30.169138888, longi=-70.805888, alt=2379.5, zmeas=0.0, SolSystemTarget='Sun')

        if np.isclose(a=result6[0], b=819.4474, atol=1e-2, rtol=0):
            a.append('result6')
            b+=1

        result7 = utc_tdb.JDUTC_to_HJDTDB(JDUTC=2458000, obsname='KPNO')


        if np.isclose(a=result7[0], b=2457999.99497543, atol=1e-7, rtol=0):
            a.append('result7')
            b+=1

        if b==8:
            print('***********SUCCESS**************\nAll barycentric correction velocities,  and time stamp conversions match expected values.\n')
        else:
            print('{} out of 8 results match. Compare outputs vs those on the github wiki. Check others - \n'.format(b,a))

        return result, result2, result3, result4, JDUTCMID, warning4, status4, corr_time, result5, result6, result7

if __name__ == "__main__":
    run_sample()
