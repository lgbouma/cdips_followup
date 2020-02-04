"""
one-time inserts to ephemerides.csv
"""
from cdips_followup.manage_ephemerides import insert_ephemeris
import os
from glob import glob
import numpy as np, pandas as pd

joel_update = 1
exofoptess_insert = 0
toi_insert = 0

##################################################
# if you have joel's updates (manual datestring) #
##################################################
if joel_update:
    ephem_files = np.sort(glob(
        '../../data/updated_ephemerides/20200204/*updateephem.txt'
    ))

    for e in ephem_files:
        insert_ephemeris(ephemsourcefile=e, ephemeris_type='hartmanupdate')


###########################################################################
# if you have a list of CDIPS candidates from exofoptess, want CTOI ephem #
###########################################################################
idpath = '../../data/updated_ephemerides/20200114_targetid_list.txt'
if exofoptess_insert:
    with open(idpath, 'r') as f:
        targetid_list = f.readlines()

    for targetid in targetid_list:
        insert_ephemeris(targetid='TIC'+targetid.replace('\n',''),
                         ephemeris_type='exofoptess_ctoi')

########################################################################
# if you have a list of TOI candidates from exofoptess, want TOI ephem #
########################################################################
idpath = '../../data/updated_ephemerides/20200115_toiid_list.txt'

if toi_insert:
    with open(idpath, 'r') as f:
        targetid_list = f.readlines()

    for targetid in targetid_list:
        insert_ephemeris(targetid=targetid.replace('\n',''),
                         ephemeris_type='exofoptess_toi')
