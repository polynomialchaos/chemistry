################################################################################
# @file parse_input.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
from .parse_mechanism import parse_mechanism
from .parse_thermos import parse_thermos
from .parse_transports import parse_transports


def parse_input(mech_path, thermo_path=None, transport_path=None):
    """Parse the provided paths."""
    tmp = parse_mechanism(mech_path)

    if thermo_path:
        tmp.thermos += parse_thermos(thermo_path)

    if transport_path:
        tmp.transports += parse_transports(transport_path)

    return tmp
