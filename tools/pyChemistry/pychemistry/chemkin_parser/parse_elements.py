################################################################################
# @file parse_elements.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import logging
from pychemistry.utilities import ElementContainer, Element
from .parse_chemkin import chemkin_format_reader


def parse_elements(path, start_keys=['ELEMENTS', 'ELEM'], end_keys=['END']):
    """Parse the elements section for a given list of strings."""
    logging.info('Parse elements from path "{:}"'.format(path))

    strings = chemkin_format_reader(
        path, start_keys=start_keys, end_keys=end_keys)
    elements = ElementContainer()

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug('Found "ALL" flag')
        elements.all_flag = True
        strings = strings[1:]

    # parse elements data
    try:
        for _, (string, _) in enumerate(strings):
            for tmp_string in string.replace('/', ' ').split():
                if tmp_string.isalpha():
                    tmp_symbol = tmp_string
                    logging.debug('Add element "{:}"'.format(tmp_symbol))
                    elements[tmp_symbol] = Element(
                        symbol=tmp_symbol
                    )
                else:
                    logging.debug(
                        'Add additional element mass "{:}"'.format(tmp_string))
                    elements[tmp_symbol].mass = tmp_string
    except:
        raise Exception('Parse error in line "{:}"'.format(string))

    return elements
