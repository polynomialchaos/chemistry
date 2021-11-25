################################################################################
# @file utilities.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
from .constants import RM


def as_int(value, round_digits=None):
    """Return a float value if possible as integer."""
    tmp = float(value) if round_digits is None else round(value, round_digits)
    return int(value) if float(value).is_integer() else tmp


def as_short(value, round_digits=None):
    """Return a string if the value is not 1."""
    tmp = as_int(value, round_digits=round_digits)
    return str(tmp) if tmp != 1 else ''


def chunk_list(elements, n_chunks):
    """Return n chunks from the given list."""
    for i in range(0, len(elements), n_chunks):
        yield elements[i:i+n_chunks]


def is_number(value):
    """Retrun True/False depending on the provided value."""
    if not type(value) in (int, float, str):
        return False

    try:
        float(value)
        return True
    except ValueError:
        return False
