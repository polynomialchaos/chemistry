####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import logging, re
from collections import OrderedDict

from pychemistry.utilities import ThermoContainer, Thermo, REF_ELEMENTS, is_number
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
def parse_thermos( path, start_keys=['THERMO', 'THER'], end_keys=['END'] ):
    """Parse the thermo section for a given list of strings."""
    logging.info( 'Parse thermo datas from path "{:}"'.format( path ) )

    strings         = chemkin_format_reader( path, start_keys=start_keys, end_keys=end_keys )
    thermos         = ThermoContainer()
    thermos.bounds  = [300.000, 1000.000, 5000.000]

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug( 'Found "ALL" flag' )
        thermos.all_flag = True
        strings = strings[1:]

    # check if temperature ranges are given or not
    try:
        tmp = [float( x ) for x in strings[0][0].split()]
        thermos.bounds = sorted( tmp )
        strings = strings[1:]
    except:
        pass

    logging.debug( 'Default temperature bounds "{:}"'.format( thermos.bounds ) )

    # parse thermo data
    try:
        parse_extra = False
        for _, (string, _) in enumerate( strings ):
            # check if tabulator character is in line
            if '\t' in string:
                raise Exception( r'Line contains tabulator (\t) character.' )

            # check if the line number is provided at the 80th position
            try:
                counter = 999 if parse_extra else int( string[79:80] )
            except:
                raise Exception( 'No line integer found (length={:})'.format( len( string.strip() ) ) )

            if string.strip()[-1] == '&': parse_extra = True

            if counter == 1:
                try:
                    tmp_symbol      = string[0:18].split( ' ' ) # remove to long comment section and append data to info
                    tmp_info        = ' '.join( tmp_symbol[1:] ).strip() + string[18:24].strip()
                    tmp_symbol      = tmp_symbol[0]

                    tmp             = string[24:44].strip()
                    tmp             = [(tmp[idx:idx+2].strip(), tmp[idx+2:idx+5].strip()) for idx in range( 0, len( tmp ), 5 )]
                    tmp_composition = OrderedDict( (key, float( nu )) for key, nu in tmp if key and float( nu ) != 0.0 )

                    tmp_phase       = string[44:45].strip() or 'G'
                    tmp_bounds      = [
                        float( string[45:55].strip().split( ' ' )[0].strip() or thermos.bounds[0] ),
                        float( string[55:65].strip().split( ' ' )[0].strip() or thermos.bounds[2] ),
                        float( string[65:78].strip().split( ' ' )[0].strip() or thermos.bounds[1] ),
                    ]

                    logging.debug( 'Values "{:}"'.format( [tmp_symbol, tmp_info, tmp_composition, tmp_phase, tmp_bounds] ) )
                except:
                    raise Exception( 'Invalid data in line counter "{:}"'.format( counter ) )
            elif counter == 999:
                # add additional element information provided after & delimiter
                parse_extra = False
                data = string.split( ' ' )
                for el, nA in zip( data[::2], data[1::2] ):
                    if el in REF_ELEMENTS and is_number( nA ):
                        tmp_composition[el] = float( nA )

                logging.debug( 'Updated elemental composition "{:}"'.format( tmp_composition ) )
            elif counter == 2:
                tmp_coeff_high  = [float( string[idx:idx+15] ) for idx in range( 0, 75, 15 )]
            elif counter == 3:
                tmp_coeff_high += [float( string[idx:idx+15] ) for idx in range( 0, 30, 15 )]
                tmp_coeff_low   = [float( string[idx:idx+15] ) for idx in range( 30, 75, 15 )]
                logging.debug( 'Coeff Low "{:}"'.format( tmp_coeff_high ) )
            elif counter == 4:
                tmp_coeff_low  += [float( string[idx:idx+15] ) for idx in range( 0, 60, 15 )]
                logging.debug( 'Coeff High "{:}"'.format( tmp_coeff_low ) )
            else:
                raise( Exception( 'Got unsupported line counter "{:}"!'.format( counter ) ) )

            if counter == 4:
                logging.debug( 'Add thermo data for species "{:}"'.format( tmp_symbol ) )

                if tmp_symbol not in thermos:
                    thermos[tmp_symbol] = Thermo(
                        symbol      = tmp_symbol,
                        info        = tmp_info,
                        composition = tmp_composition,
                        phase       = tmp_phase,
                        bounds      = tmp_bounds,
                        coeff_low   = tmp_coeff_low,
                        coeff_high  = tmp_coeff_high
                    )
    except:
        logging.error( 'Parse error in line "{:}"'.format( string ) )
        raise

    return thermos