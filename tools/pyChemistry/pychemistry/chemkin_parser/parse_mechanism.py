####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import logging, os

from pychemistry.utilities import Mechanism
from .parse_elements import parse_elements
from .parse_specii import parse_specii
from .parse_reactions import parse_reactions
from .parse_thermos import parse_thermos
from .parse_transports import parse_transports

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
def parse_mechanism( path ):
    """Parse the provided path."""
    logging.info( 'Parse mechanism from path "{:}"'.format( path ) )

    return Mechanism(
        name        = os.path.basename( path ),
        elements    = parse_elements( path ),
        specii      = parse_specii( path ),
        reactions   = parse_reactions( path ),
        thermos     = parse_thermos( path ),
        transports  = parse_transports( path, start_keys=['TRANSPORT', 'TRANS'], end_keys=['END'] )
    )