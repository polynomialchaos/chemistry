################################################################################
# @file parse_reactions.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import logging
import re
from pychemistry.utilities import ReactionContainer, Reaction, ReactionType
from pychemistry.utilities import conv_k0, conv_Ea
from pychemistry.utilities import cmsk_to_si
from .parse_chemkin import chemkin_format_reader

regex_nu = re.compile(r'^([-+]?((\d+\.\d*)|(\.\d+)|(\d+))([eE][-+]?\d+)?)')
regex_pressure = re.compile(r'\(\+.*?\)')
regex_three = re.compile(r'\+\s*M')
ndef_keys = ['SRI', 'LT', 'JAN', 'FIT1', 'HV',
             'TDEP', 'EXCI', 'MOME', 'XSMI', 'PLOG', 'UNITS']


def parse_auxiliary_data(strings):
    """Parse auxiliary data for a given list of strings."""
    result = {}

    for _, string in enumerate(strings):
        # split the provided string into its parts
        tmpString = [y.strip() for x in string.split('/')
                     for y in x.split() if len(y) != 0]

        # if only one value add 'FLAG' key
        if len(tmpString) == 1:
            tmpString = ['FLAGS'] + tmpString

        # search for unsupported keys
        key = tmpString[0]
        if key in ndef_keys:
            raise KeyError(
                'Unsupported auxiliary key "{:}" provided!'.format(key))

        # depending on the keyword store the data in the auxiliary dict
        if key in ['LOW', 'HIGH', 'TROE', 'REV']:
            if key in result:
                raise KeyError(
                    'Redefinition of auxiliary key "{:}"!'.format(key))

            result[key] = [float(x) for x in tmpString[1:]]
        elif key in ['FORD', 'RORD']:
            if not key in result:
                result[key] = {}

            for i in range(1, len(tmpString[1:]), 2):
                result[key][tmpString[i]] = float(tmpString[i+1])
        elif key in ['FLAGS']:
            if not key in result:
                result[key] = []

            result[key].append(tmpString[1])
        elif key in ['UNITS']:
            if not key in result:
                result[key] = []

            result[key].extend(tmpString[1:])
        else:
            key = 'EFFS'
            if not key in result:
                result[key] = {}

            for i in range(0, len(tmpString), 2):
                result[key][tmpString[i]] = float(tmpString[i+1])

        logging.debug('Auxiliary data "%s" = "%s"', key, result[key])

    return result


def parse_reaction_line(string):
    """Parse reactants and products strings."""
    # remove unwanted characters and split into parts
    tmp_string = regex_pressure.sub('', string)
    tmp_string = regex_three.sub('', tmp_string)
    tmp_string = tmp_string.split('+')

    # split into nu and species
    # (take care of multiple species definitions)
    result = {}
    for tmp in tmp_string:
        match = regex_nu.search(tmp)
        idx = 0 if match is None else match.span()[-1]
        match_nu = float(tmp[:idx]) if tmp[:idx] else 1.0
        match_sp = tmp[idx:].strip()

        if match_sp in result:
            result[match_sp] += match_nu
        else:
            result[match_sp] = match_nu

    return result


def parse_reaction(strings, unit_k0, unit_ea):
    """Parse a reaction for a given list of strings."""
    tmp_reaction_line = strings[0]
    tmp_auxiliary_lines = strings[1:]
    logging.debug('Add reaction "%s"', tmp_reaction_line)
    for string in tmp_auxiliary_lines:
        logging.debug('Additional line "%s"', string)

    # reaction type and falloff species
    if regex_pressure.search(tmp_reaction_line):
        tmp_type = ReactionType.PRESSURE
        tmp_falloff_species = tmp_reaction_line.split(
            '(+')[1].split(')')[0].strip()
    elif regex_three.search(tmp_reaction_line):
        tmp_type = ReactionType.THREE_BODY
        tmp_falloff_species = 'M'
    else:
        tmp_type = ReactionType.DEFAULT
        tmp_falloff_species = ''

    logging.debug('Reaction type: %s', tmp_type)

    # arrhenius coefficients
    tmp_arr_coeff = [float(x) for x in tmp_reaction_line.split()[-3:]]
    logging.debug('Arrhenius coefficients: "%s"', tmp_arr_coeff)

    # reversibility
    tmp_line = ''.join(tmp_reaction_line.split()[:-3])
    if ('=>' in tmp_line) and (not '<=>' in tmp_line):
        tmp_delimiter = '=>'
        tmp_is_reversible = False
    elif '<=>' in tmp_line:
        tmp_delimiter = '<=>'
        tmp_is_reversible = True
    elif '=' in tmp_line:
        tmp_delimiter = '='
        tmp_is_reversible = True
    else:
        raise Exception(
            'No supported delimiter string "{:}" provided!'.format(tmp_line))

    logging.debug('Reversibility: %s', tmp_is_reversible)

    # reactants and products
    tmp_reactants = parse_reaction_line(tmp_line.split(tmp_delimiter)[0])
    tmp_products = parse_reaction_line(tmp_line.split(tmp_delimiter)[1])

    logging.debug('Reactants "%s"', tmp_reactants)
    logging.debug('Products "%s"', tmp_products)

    # auxiliary data
    tmp_auxiliary = parse_auxiliary_data(tmp_auxiliary_lines)
    tmp_adv_arr_key = 'HIGH' if 'HIGH' in tmp_auxiliary else 'LOW'

    reaction = Reaction(reaction_type=tmp_type, reactants=tmp_reactants,
                        products=tmp_products, is_reversible=tmp_is_reversible)

    reaction.falloff_species = tmp_falloff_species
    reaction.adv_arr_key = tmp_adv_arr_key

    if 'TROE' in tmp_auxiliary:
        reaction.troe_coeff = tmp_auxiliary['TROE']

    if 'FORD' in tmp_auxiliary:
        reaction.f_orders = tmp_auxiliary['FORD']

    if 'RORD' in tmp_auxiliary:
        reaction.r_orders = tmp_auxiliary['RORD']

    if 'FLAGS' in tmp_auxiliary:
        reaction.flags = tmp_auxiliary['FLAGS']

    if 'EFFS' in tmp_auxiliary:
        reaction.efficiencies = tmp_auxiliary['EFFS']

    reaction.arr_coeff = cmsk_to_si(
        tmp_arr_coeff, reaction.f_conv_si(), unit_k0, unit_ea)

    if 'REV' in tmp_auxiliary:
        reaction.rev_arr_coeff = cmsk_to_si(
            tmp_auxiliary['REV'], reaction.r_conv_si(), unit_k0, unit_ea)

    if 'LOW' in tmp_auxiliary:
        reaction.adv_arr_coeff = cmsk_to_si(
            tmp_auxiliary[tmp_adv_arr_key],
            reaction.f_conv_si(add=1.0),
            unit_k0, unit_ea
        )
    elif 'HIGH' in tmp_auxiliary:
        raise NotImplementedError(
            'HIGH keyword is not supported. ' +
            'If use, check add parameter in unit conversion (may be 0)')
        # reaction.adv_arr_coeff = cmsk_to_si(
        #     tmp_auxiliary[tmp_adv_arr_key],
        #     reaction.f_conv_si(add=-1.0), unit_k0, unit_ea)

    return reaction


def parse_reactions(path, start_keys=None, end_keys=None):
    """Parse the reactions section for a given list of strings."""
    logging.info('Parse reactions from path "%s"', path)
    start_keys = ['REACTIONS', 'REAC'] if start_keys is None else start_keys
    end_keys = ['END'] if end_keys is None else end_keys

    strings = chemkin_format_reader(
        path, start_keys=start_keys, end_keys=end_keys)
    reactions = ReactionContainer()

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug('Found "ALL" flag')
        reactions.all_flag = True
        strings = strings[1:]

    # check if units are provided
    if not '=' in strings[0][0]:
        for key in strings[0][0].split():
            if key in conv_k0.keys():
                reactions.unit_k0 = key
            elif key in conv_Ea.keys():
                reactions.unit_ea = key
            else:
                raise KeyError(
                    'Provided unsupported reaction unit: "{:}"!'.format(key))

        strings = strings[1:]

    logging.debug('Reaction default units "%s" and "%s"',
                  reactions.unit_k0, reactions.unit_ea)

    # parse reactions data
    try:
        tmp_reaction = []
        for _, (string, _) in enumerate(strings):
            # check if reaction tmp_delimiter is found
            if '=' in string and tmp_reaction:
                reactions.append(
                    parse_reaction(
                        tmp_reaction, reactions.unit_k0, reactions.unit_ea)
                )

                # free the tmpReactions
                tmp_reaction = []

            # append line if no special character is found
            tmp_reaction.append(string)

        # append last reaction
        if tmp_reaction:
            reactions.append(
                parse_reaction(tmp_reaction, reactions.unit_k0,
                               reactions.unit_ea)
            )
    except:
        logging.error('Parse error in line "%s"', string)
        raise

    return reactions
