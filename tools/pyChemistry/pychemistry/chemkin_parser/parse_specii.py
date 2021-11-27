################################################################################
# @file parse_specii.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import logging
from pychemistry.utilities import SpeciesContainer, Species
from .parse_chemkin import chemkin_format_reader


def parse_specii(path, start_keys=None, end_keys=None):
    """Parse the species section for a given list of strings."""
    logging.info('Parse species from path "%s"', path)
    start_keys = ['SPECIES', 'SPEC'] if start_keys is None else start_keys
    end_keys = ['END'] if end_keys is None else end_keys

    strings = chemkin_format_reader(
        path, start_keys=start_keys, end_keys=end_keys)
    species = SpeciesContainer()

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug('Found "ALL" flag')
        species.all_flag = True
        strings = strings[1:]

    # parse species data
    try:
        for _, (string, _) in enumerate(strings):
            for tmp_string in string.split():
                logging.debug('Add species "%s"', tmp_string)

                species[tmp_string] = Species(
                    symbol=tmp_string
                )
    except:
        logging.error('Parse error in line "%s"', string)
        raise

    return species
