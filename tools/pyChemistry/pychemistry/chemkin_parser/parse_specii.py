####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import logging
from pychemistry.utilities import SpeciesContainer, Species
from .parse_chemkin import chemkin_format_reader

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
def parse_specii( path, start_keys=['SPECIES', 'SPEC'], end_keys=['END'] ):
    """Parse the species section for a given list of strings."""
    logging.info( 'Parse species from path "{:}"'.format( path ) )

    strings = chemkin_format_reader( path, start_keys=start_keys, end_keys=end_keys )
    species = SpeciesContainer()

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug( 'Found "ALL" flag' )
        species.all_flag = True
        strings = strings[1:]

    # parse species data
    try:
        for _, (string, _) in enumerate( strings ):
            for tmp_string in string.split():
                logging.debug( 'Add species "{:}"'.format( tmp_string ) )

                species[tmp_string] = Species(
                    symbol  = tmp_string
                )
    except:
        logging.error( 'Parse error in line "{:}"'.format( string ) )
        raise

    return species