import pandas as pd

from cdips.utils.catalogs import get_exofop_toi_catalog

def ticid_to_toiid(tic_id):

    assert isinstance(tic_id, str)

    toidf = get_exofop_toi_catalog()

    sel = toidf['TIC ID'].astype(str) == tic_id

    try:
        toi_id = str(toidf[sel]['TOI'].iloc[0])
        if len(toi_id) > 1:
            return str(toi_id)
        else:
            return None
    except:
        return None
