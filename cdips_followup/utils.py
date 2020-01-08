from cdips.utils.catalogs import ticid_to_toiid as cdips_ticid_to_toiid

def ticid_to_toiid(tic_id):

    return cdips_ticid_to_toiid(tic_id)


def ticid_and_toiid_to_targetid(tic_id, toi_id):
    """
    If TOI identifier exists, use that. Otherwise, try for TIC identifer. If
    neither was passed, return None.
    """

    if isinstance(tic_id, str) and isinstance(toi_id, str):

        return toi_id + '.01'

    elif isinstance(tic_id, str) and toi_id is None:

        return tic_id + '.01'

    else:

        return None
