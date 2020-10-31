####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import numpy as np
from collections import OrderedDict

from .constants import RM

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
def as_int( value, round_digits=None ):
    """Return a float value if possible as integer."""
    tmpValue = float( value ) if round_digits is None else round( value, round_digits )
    return int( value ) if float( value ).is_integer() else tmpValue

def as_short( value, round_digits=None ):
    """Return a string if the value is not 1."""
    return str( as_int( value, round_digits=round_digits ) ) if value != 1 else ''

def chunk_list( elements, n_chunks ):
    """Return n chunks from the given list."""
    for i in range( 0, len( elements ), n_chunks ):
        yield elements[i:i+n_chunks]

def f_arr_ta( T, coeff ):
    """Return the evaluation of the arrhenius function (using activation temperature)."""
    return f_arr( T, coeff, conv=1.0 )

def f_arr_ea( T, coeff ):
    """Return the evaluation of the arrhenius function (using activation energy)."""
    return f_arr( T, coeff, conv=RM )

def f_arr( T, coeff, conv=1.0 ):
    """Return the evaluation of the arrhenius function (using conversion)."""
    return coeff[0] * np.exp( coeff[1] * np.log( T ) - coeff[2] / ( T * conv ) )

def alpha_to_arr( c_alpha, T_alpha_max, n_alpha ):
    """Return the arrhenius form of the alpha coefficients."""
    return [c_alpha * T_alpha_max**n_alpha * np.exp( n_alpha ), -n_alpha, n_alpha * T_alpha_max]

def is_number( value ):
    """Retrun True/False depending on the provided value."""
    try:
        float( value )
        return True
    except ValueError:
        return False