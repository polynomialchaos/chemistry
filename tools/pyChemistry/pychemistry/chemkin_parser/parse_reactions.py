####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import logging, re
from collections import OrderedDict

from pychemistry.utilities import ReactionContainer, Reaction, ReactionType, conv_k0, def_unit_k0, conv_Ea, def_unit_Ea, cmsk_to_si
from .parse_chemkin import chemkin_format_reader

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
regex_pressure  = re.compile( r'\(\+.*?\)' )
regex_three     = re.compile( r'\+\s*M' )
ndef_keys       = ['SRI', 'LT', 'JAN', 'FIT1', 'HV', 'TDEP', 'EXCI', 'MOME', 'XSMI', 'PLOG', 'UNITS']

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
def parse_auxiliary_data( strings ):
    """Parse auxiliary data for a given list of strings."""
    result = OrderedDict()

    for idx, string in enumerate( strings ):
        # split the provided string into its parts
        tmpString = [y.strip() for x in string.split( '/' ) for y in x.split() if len( y ) != 0]

        # if only one value add 'FLAG' key
        if len( tmpString ) == 1:
            tmpString = ['FLAGS'] + tmpString

        # search for unsupported keys
        key = tmpString[0]
        if key in ndef_keys:
            raise( KeyError( 'Unsupported auxiliary key "{:}" provided!'.format( key ) ) )

        # depending on the keyword store the data in the auxiliary dict
        if key in ['LOW', 'HIGH', 'TROE', 'REV']:
            if key in result:
                raise( KeyError( 'Redefinition of auxiliary key "{:}"!'.format( key ) ) )

            result[key] = [float( x ) for x in tmpString[1:]]
        elif key in ['FORD', 'RORD']:
            if not key in result:
                result[key] = OrderedDict()

            for idx in range( 1, len( tmpString[1:] ), 2 ):
                result[key][tmpString[idx]] = float( tmpString[idx+1] )
        elif key in ['FLAGS']:
            if not key in result:
                result[key] = []

            result[key].append( tmpString[1] )
        elif key in ['UNITS']:
            if not key in result:
                result[key] = []

            result[key].extend( tmpString[1:] )
        else:
            key = 'EFFS'
            if not key in result:
                result[key] = OrderedDict()

            for idx in range( 0, len( tmpString ), 2 ):
                result[key][tmpString[idx]] = float( tmpString[idx+1] )

        logging.debug( 'Auxiliary data "{:}" = "{:}"'.format( key, result[key] ) )

    return result

def parse_reaction_line( string ):
    """Parse reactants and products strings."""

    # remove unwanted characters and split into parts
    tmp_string = regex_pressure.sub( '', string )
    tmp_string = regex_three.sub( '', tmp_string )
    tmp_string = tmp_string.split( '+' )

    # split into nu and species
    # (take care of multiple species definitions)
    result = OrderedDict()
    for tmp in tmp_string:
        match       = re.search( r'[^\W\d]', tmp )
        match_sp    = tmp[match.start():].strip()
        match_nu    = float( tmp[:match.start()] ) if tmp[:match.start()] else 1.0

        if match_sp in result:
            result[match_sp] += match_nu
        else:
            result[match_sp] = match_nu

    return result

def parse_reaction( strings, unit_k0, unit_Ea ):
    """Parse a reaction for a given list of strings."""
    tmp_reaction_line = strings[0]
    logging.debug( 'Add reaction "{:}"'.format( tmp_reaction_line ) )

    tmp_auxiliary_lines = strings[1:]
    for string in tmp_auxiliary_lines:
        logging.debug( 'Additional line "{:}"'.format( string ) )

    # reaction type and falloff species
    if regex_pressure.search( strings[0] ):
        tmp_type            = ReactionType.PRESSURE
        tmp_falloff_species = strings[0].split( '(+' )[1].split( ')' )[0].strip()
    elif regex_three.search( strings[0] ):
        tmp_type            = ReactionType.THREE_BODY
        tmp_falloff_species = 'M'
    else:
        tmp_type            = ReactionType.DEFAULT
        tmp_falloff_species = ''

    logging.debug( 'Reaction type: {:}'.format( tmp_type ) )

    # arrhenius coefficients
    tmp_arr_coeff = [float( x ) for x in tmp_reaction_line.split()[-3:]]
    logging.debug( 'Arrhenius coefficients: {:}'.format( tmp_arr_coeff ) )

    # reversibility
    tmp_line = ''.join( tmp_reaction_line.split()[:-3] )
    if ('=>' in tmp_line) and (not '<=>' in tmp_line):
        tmp_delimiter       = '=>'
        tmp_is_reversible   = False
    elif '<=>' in tmp_line:
        tmp_delimiter       = '<=>'
        tmp_is_reversible   = True
    elif '=' in tmp_line:
        tmp_delimiter       = '='
        tmp_is_reversible   = True
    else:
        raise( Exception( 'No supported delimiter string "{:}" provided!'.format( tmp_line ) ) )

    logging.debug( 'Reversibility: {:}'.format( tmp_is_reversible ) )

    # reactants and products
    tmp_reactants   = parse_reaction_line( tmp_line.split( tmp_delimiter )[0] )
    tmp_products    = parse_reaction_line( tmp_line.split( tmp_delimiter )[1] )

    logging.debug( 'Reactants "{:}"'.format( tmp_reactants ) )
    logging.debug( 'Products "{:}"'.format( tmp_products ) )

    # auxiliary data
    tmp_auxiliary   = parse_auxiliary_data( tmp_auxiliary_lines )
    tmp_adv_arr_key = 'HIGH' if 'HIGH' in tmp_auxiliary else 'LOW'

    reaction = Reaction(
        reaction_type   = tmp_type,
        reactants       = tmp_reactants,
        products        = tmp_products,
        is_reversible   = tmp_is_reversible,
        falloff_species = tmp_falloff_species,
        adv_arr_key     = tmp_adv_arr_key,
        troe_coeff      = tmp_auxiliary.get( 'TROE', None ),
        f_orders        = tmp_auxiliary.get( 'FORD', None ),
        r_orders        = tmp_auxiliary.get( 'RORD', None ),
        flags           = tmp_auxiliary.get( 'FLAGS', None ),
        efficiencies    = tmp_auxiliary.get( 'EFFS', None )
    )

    reaction.arr_coeff = cmsk_to_si( tmp_arr_coeff, reaction.f_conv_si(), unit_k0, unit_Ea )

    if 'REV' in tmp_auxiliary:
        reaction.rev_arr_coeff = cmsk_to_si( tmp_auxiliary['REV'], reaction.r_conv_si(), unit_k0, unit_Ea )

    if 'LOW' in tmp_auxiliary:
        reaction.adv_arr_coeff = cmsk_to_si( tmp_auxiliary[tmp_adv_arr_key], reaction.f_conv_si( add=1.0 ), unit_k0, unit_Ea )
    elif 'HIGH' in tmp_auxiliary:
        reaction.adv_arr_coeff = cmsk_to_si( tmp_auxiliary[tmp_adv_arr_key], reaction.f_conv_si( add=-1.0 ), unit_k0, unit_Ea )

    return reaction

def parse_reactions( path, start_keys=['REACTIONS', 'REAC'], end_keys=['END'] ):
    """Parse the reactions section for a given list of strings."""
    logging.info( 'Parse reactions from path "{:}"'.format( path ) )

    strings     = chemkin_format_reader( path, start_keys=start_keys, end_keys=end_keys )
    reactions   = ReactionContainer()

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug( 'Found "ALL" flag' )
        reactions.all_flag = True
        strings = strings[1:]

    # check if units are provided
    if not '=' in strings[0][0]:
        for key in [x.strip() for x in strings[0][0].split( ' ' ) if x]:
            if key in conv_k0.keys():
                reactions.unit_k0 = key
            elif key in conv_Ea.keys():
                reactions.unit_Ea = key
            else:
                raise( KeyError( 'Provided unsupported reaction unit: "{:}"!'.format( key ) ) )

        strings = strings[1:]

    logging.debug( 'Reaction default units "{:}" and "{:}"'.format( reactions.unit_k0, reactions.unit_Ea ) )

    # parse reactions data
    try:
        tmp_reaction = []
        for _, (string, _) in enumerate( strings ):
            # check if reaction tmp_delimiter is found
            if '=' in string and tmp_reaction:
                reactions.append(
                    parse_reaction( tmp_reaction, reactions.unit_k0, reactions.unit_Ea )
                )

                # free the tmpReactions
                tmp_reaction = []

            # append line if no special character is found
            tmp_reaction.append( string )

        # append last reaction
        if tmp_reaction:
            reactions.append(
                parse_reaction( tmp_reaction, reactions.unit_k0, reactions.unit_Ea )
            )
    except:
        logging.error( 'Parse error in line "{:}"'.format( string ) )
        raise

    return reactions
