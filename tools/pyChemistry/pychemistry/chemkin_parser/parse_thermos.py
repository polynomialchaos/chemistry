################################################################################
# @file parse_thermos.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import logging
from pychemistry.utilities import ThermoContainer, Thermo
from pychemistry.utilities import REF_ELEMENTS, is_number
from .parse_chemkin import chemkin_format_reader


def parse_thermos(path, start_keys=None, end_keys=None):
    """Parse the thermo section for a given list of strings."""
    logging.info('Parse thermo datas from path "%s"', path)
    start_keys = ['THERMO', 'THER'] if start_keys is None else start_keys
    end_keys = ['END'] if end_keys is None else end_keys

    strings = chemkin_format_reader(
        path, start_keys=start_keys, end_keys=end_keys)
    thermos = ThermoContainer()
    thermos.bounds = [300.000, 1000.000, 5000.000]

    # check if ALL flag is provided
    if strings and 'ALL' in strings[0][0]:
        logging.debug('Found "ALL" flag')
        thermos.all_flag = True
        strings = strings[1:]

    # check if temperature ranges are given or not
    try:
        tmp = [float(x) for x in strings[0][0].split()]
        thermos.bounds = sorted(tmp)
        strings = strings[1:]
    except IndexError:
        pass

    logging.debug('Default temperature bounds "%s"', thermos.bounds)

    # parse thermo data
    try:
        parse_extra = False
        for _, (string, _) in enumerate(strings):
            # check if tabulator character is in line
            if '\t' in string:
                raise Exception('Line contains tabulator (\\t) character.')

            # check if the line number is provided at the 80th position
            try:
                counter = 999 if parse_extra else int(string[79:80])
            except ValueError as missing_index:
                raise Exception(
                    'No line integer found (length={:})'.format(
                        len(string.strip()))) from missing_index

            if string.strip()[-1] == '&':
                parse_extra = True

            if counter == 1:
                # remove to long comment section and append data to info
                tmp_symbol = string[0:18].split()
                tmp_info = ' '.join(
                    tmp_symbol[1:]).strip() + string[18:24].strip()
                tmp_symbol = tmp_symbol[0]

                tmp = string[24:44].strip()
                tmp = [(tmp[idx:idx+2].strip(), tmp[idx+2:idx+5].strip())
                       for idx in range(0, len(tmp), 5)]
                tmp_composition = {key: float(
                    value) for key, value in tmp if key and float(value) != 0}

                tmp_phase = string[44:45].strip() or 'G'
                tmp_bounds = [
                    float(string[45:55].split()[0] or thermos.bounds[0]),
                    float(string[55:65].split()[0] or thermos.bounds[2]),
                    float(string[65:78].split()[0] or thermos.bounds[1]),
                ]

                logging.debug('Values "%s"', [tmp_symbol, tmp_info,
                                              tmp_composition, tmp_phase,
                                              tmp_bounds])
            elif counter == 999:
                # add additional element information provided after & delimiter
                parse_extra = False
                data = string.split()
                if (len(data) % 2) != 0:
                    raise Exception('Invalid number of additional '
                                    'composition strings!')

                for ael, ael_n in zip(data[::2], data[1::2]):
                    if ael in REF_ELEMENTS and is_number(ael_n):
                        tmp_composition[ael] = float(ael_n)
                    else:
                        raise Exception('Invalid element data in line!')

                logging.debug(
                    'Updated elemental composition "%s"', tmp_composition)
            elif counter == 2:
                tmp_coeff_high = [float(string[idx:idx+15])
                                  for idx in range(0, 75, 15)]
            elif counter == 3:
                tmp_coeff_high += [float(string[idx:idx+15])
                                   for idx in range(0, 30, 15)]
                tmp_coeff_low = [float(string[idx:idx+15])
                                 for idx in range(30, 75, 15)]
                logging.debug('Coeff Low "%s"', tmp_coeff_high)
            elif counter == 4:
                tmp_coeff_low += [float(string[idx:idx+15])
                                  for idx in range(0, 60, 15)]
                logging.debug('Coeff High "%s"', tmp_coeff_low)
            else:
                raise Exception('Got unsupported line counter "{:}"!'.format(
                    counter))

            if counter == 4:
                logging.debug(
                    'Add thermo data for species "%s"', tmp_symbol)

                if tmp_symbol not in thermos:
                    thermos[tmp_symbol] = Thermo(tmp_symbol)
                    thermos[tmp_symbol].info = tmp_info
                    thermos[tmp_symbol].composition = tmp_composition
                    thermos[tmp_symbol].phase = tmp_phase
                    thermos[tmp_symbol].bounds = tmp_bounds
                    thermos[tmp_symbol].coeff_low = tmp_coeff_low
                    thermos[tmp_symbol].coeff_high = tmp_coeff_high
                else:
                    logging.warning(
                        'Ignore redefinition of thermo data for "%s"!',
                        tmp_symbol)

    except:
        logging.error('Parse error in line "%s"', string)
        raise

    return thermos
