################################################################################
# @file parse_mechanism.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import os
import logging
from pychemistry.utilities import Mechanism
from .parse_elements import parse_elements
from .parse_specii import parse_specii
from .parse_reactions import parse_reactions
from .parse_thermos import parse_thermos
from .parse_transports import parse_transports


def parse_mechanism(path):
    """Parse the provided path."""
    logging.info('Parse mechanism from path "%s"', path)

    return Mechanism(
        name=os.path.basename(path),
        elements=parse_elements(path),
        specii=parse_specii(path),
        reactions=parse_reactions(path),
        thermos=parse_thermos(path),
        transports=parse_transports(
            path, start_keys=['TRANSPORT', 'TRANS'], end_keys=['END'])
    )
