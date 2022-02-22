################################################################################
# @file parse_transports.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import logging
from pychemistry.utilities import TransportContainer, Transport, CGS_TO_SI_COLLJ
from pychemistry.utilities import CGS_TO_SI_DIPMO, CGS_TO_SI_POL
from .parse_chemkin import chemkin_format_reader


def parse_transports(path, start_keys=None, end_keys=None):
    """Parse the transport section for a given list of strings."""
    logging.info('Parse transport datas from path "%s"', path)

    strings = chemkin_format_reader(
        path, start_keys=start_keys, end_keys=end_keys)
    transports = TransportContainer()

    # check if TRANS flag is provided
    if strings:
        if not start_keys and 'TRANS' in strings[0][0]:
            strings = strings[1:]
            logging.debug('Removed "TRANS" flag in transport file')

    # check if END flag is provided
    if strings:
        if not end_keys and 'END' in strings[-1][0]:
            strings = strings[:-1]
            logging.debug('Removed "END" flag in transport file')

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug('Found "ALL" flag')
        transports.all_flag = True
        strings = strings[1:]

    # parse transport data
    try:
        for _, (string, _) in enumerate(strings):
            tmp_string = string.split()
            logging.debug(
                'Add transport data for species "%s"', tmp_string[0])
            logging.debug('Values "%s"', tmp_string[1:])

            if tmp_string[0] not in transports:
                transports[tmp_string[0]] = Transport(tmp_string[0])
                transports[tmp_string[0]].geom = int(tmp_string[1])
                transports[tmp_string[0]].pot_lj = float(tmp_string[2])
                transports[tmp_string[0]].col_lj = float(
                    tmp_string[3]) * CGS_TO_SI_COLLJ
                transports[tmp_string[0]].dip_mo = float(
                    tmp_string[4]) * CGS_TO_SI_DIPMO
                transports[tmp_string[0]].pol = float(
                    tmp_string[5]) * CGS_TO_SI_POL
                transports[tmp_string[0]].rot_rel = tmp_string[6]
            else:
                logging.warning(
                    'Ignore redefinition of transport data for "%s"!',
                    tmp_string[0])

    except:
        logging.error('Parse error in line "%s"', string)
        raise

    return transports
