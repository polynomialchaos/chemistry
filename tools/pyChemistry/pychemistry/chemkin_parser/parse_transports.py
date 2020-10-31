####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import logging
from pychemistry.utilities import TransportContainer, Transport, CGS_TO_SI_COLLJ, CGS_TO_SI_DIPMO, CGS_TO_SI_POL
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
def parse_transports( path, start_keys=None, end_keys=None ):
    """Parse the transport section for a given list of strings."""
    logging.info( 'Parse transport datas from path "{:}"'.format( path ) )

    strings     = chemkin_format_reader( path, start_keys=start_keys, end_keys=end_keys )
    transports  = TransportContainer()

    # check if END flag is provided
    if strings:
        if not start_keys and 'TRANS' in strings[0][0]:
            strings = strings[1:]
            logging.debug( 'Removed "TRANS" flag in transport file' )

    # check if END flag is provided
    if strings:
        if not end_keys and 'END' in strings[-1][0]:
            strings = strings[:-1]
            logging.debug( 'Removed "END" flag in transport file' )

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug( 'Found "ALL" flag' )
        transports.all_flag = True
        strings = strings[1:]

    # parse transport data
    try:
        for _, (string, _) in enumerate( strings ):
            tmp_string = string.split()
            logging.debug( 'Add transport data for species "{:}"'.format( tmp_string[0] ) )
            logging.debug( 'Values "{:}"'.format( tmp_string[1:] ) )

            if tmp_string[0] not in transports:
                transports[tmp_string[0]] = Transport(
                    symbol  = tmp_string[0],
                    geom    = int( tmp_string[1] ),
                    pot_lj  = float( tmp_string[2] ),
                    col_lj  = float( tmp_string[3] ) * CGS_TO_SI_COLLJ,
                    dip_mo  = float( tmp_string[4] ) * CGS_TO_SI_DIPMO,
                    pol     = float( tmp_string[5] ) * CGS_TO_SI_POL,
                    rot_rel = tmp_string[6],
                )
    except:
        logging.error( 'Parse error in line "{:}"'.format( string ) )
        raise

    return transports