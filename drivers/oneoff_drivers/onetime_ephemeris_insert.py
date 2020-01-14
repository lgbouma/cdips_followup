"""
one-time inserts to ephemerides.csv
"""
from cdips_followup.manage_ephemerides import insert_ephemeris
import os
from glob import glob
import numpy as np, pandas as pd

# then include joel's updates
ephem_files = np.sort(glob(
    '../../data/updated_ephemerides/20200114/*updateephem.txt'
))

for e in ephem_files:
    insert_ephemeris(ephemsourcefile=e, ephemeris_type='hartmanupdate')
