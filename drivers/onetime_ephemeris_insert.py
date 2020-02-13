"""
one-time inserts to ephemerides.csv
"""
from cdips_followup.manage_ephemerides import insert_ephemeris
import os
from glob import glob
import numpy as np, pandas as pd

joel_update = 1
exofoptess_ctoi_insert = 0
toi_insert = 0
pipe_insert = 0

##################################################
# if you have joel's updates (manual datestring) #
##################################################
if joel_update:
    ephem_files = np.sort(glob(
        '../../data/updated_ephemerides/20200211/*updateephem.txt'
    ))

    for e in ephem_files:
        insert_ephemeris(ephemsourcefile=e, ephemeris_type='hartmanupdate')


###########################################################################
# if you have a list of CDIPS candidates from exofoptess, want CTOI ephem #
###########################################################################
idpath = '../../data/updated_ephemerides/20200207_targetid_list.txt'
if exofoptess_ctoi_insert:
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

############################################################################
# if you have the output from cdipspipeline that went to exofop (e.g., for #
# updated ephemerides, rather than those initially in the CTOI table)      #
############################################################################
ephemsourcefile = '../../data/updated_ephemerides/20200207_sectors_8_thru_11_clear_threshold_w_sourceid.csv'

if pipe_insert:

    df = pd.read_csv(ephemsourcefile, sep="|")

    for targetid in np.array(df.target):

        insert_ephemeris(targetid=targetid, ephemeris_type='cdipspipeline',
                         ephemsourcefile=ephemsourcefile)