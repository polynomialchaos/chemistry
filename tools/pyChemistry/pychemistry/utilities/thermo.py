################################################################################
# @file thermo.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
from .base import Base, BaseDictContainer
from .constants import REF_ELEMENTS


PHASE_NUM = {'G': 0}    # Conversion: Phase string
N_BOUNDS = 3            # Number of temperature bounds
N_NASA = 7              # Number of NASA polynomial coefficients


class Thermo(Base):
    """Thermo object (storing data)."""

    def __init__(self, symbol, info=None, composition=None, phase=None,
                 bounds=None, coeff_low=None, coeff_high=None):
        self.symbol = symbol
        self.info = info
        self.composition = composition
        self.phase = phase
        self.bounds = bounds
        self.coeff_low = coeff_low
        self.coeff_high = coeff_high

    def __str__(self):
        return self.symbol

    def _check_list(self):
        return [
            (bool(self.composition), 'Elemental composition empty!'),
            (len(self.bounds) == N_BOUNDS,
             'NASA temp. bounds of wrong size ({:})!'.format(
                 len(self.bounds))),
            (len(self.coeff_low) == N_NASA,
             'NASA low temp. coefficients of wrong size ({:})!'.format(
                 len(self.coeff_low))),
            (len(self.coeff_high) == N_NASA,
             'NASA high temp. coefficients of wrong size ({:})!'.format(
                 len(self.coeff_high))),
            (self.phase == 0, 'Phase not valid!'),
        ]

    def _get_comp_strings(self):
        all_comps = ['{:2}{:3}'.format(
            key, int(self.composition[key])) for key in self.composition]

        if len(all_comps) > 4 or any(len(x) > 5 for x in all_comps):
            return '', ' '.join(all_comps)
        else:
            return ''.join(all_comps), ''

    @property
    def bounds(self):
        return list(sorted(self._bounds))

    @bounds.setter
    def bounds(self, value):
        if value is None:
            return
        tmp = sorted([float(x) for x in value])
        self._bounds = [tmp[0], tmp[2], tmp[1]]

    def chemkinify(self):
        symb_info = '{:} {:}'.format(self.symbol, self.info)
        short_comps, long_comps = self._get_comp_strings()
        line1 = '{:24}{:20}{:1}{:-10.3f}{:-10.3f}{:-10.3f}'.format(
            symb_info[:24], short_comps, self._phase, *self._bounds)
        line1_long = '' if not long_comps else ' &\n{:}'.format(long_comps)
        line2 = ''.join(['{: 15.8E}'.format(coeff)
                        for coeff in self.coeff_high[:5]])
        line3 = ''.join(['{: 15.8E}'.format(coeff)
                        for coeff in self.coeff_high[5:] + self.coeff_low[:3]])
        line4 = ''.join(['{: 15.8E}'.format(coeff)
                        for coeff in self.coeff_low[3:]])

        return '{:79}1{:}\n{:79}2\n{:79}3\n{:79}4'.format(
            line1, line1_long, line2, line3, line4
        )

    @property
    def coeff_high(self):
        return self._coeff_high

    @coeff_high.setter
    def coeff_high(self, value):
        if value is None:
            return
        self._coeff_high = [float(x) for x in value]

    @property
    def coeff_low(self):
        return self._coeff_low

    @coeff_low.setter
    def coeff_low(self, value):
        if value is None:
            return
        self._coeff_low = [float(x) for x in value]

    @property
    def composition(self):
        return self._composition

    @composition.setter
    def composition(self, value):
        if value is None:
            return
        self._composition = {key: float(value) for key, value in value.items()}

    @property
    def info(self):
        return self._info

    @info.setter
    def info(self, value):
        if value is None:
            return
        self._info = value.strip()

    @property
    def molar_mass(self):
        return sum(
            REF_ELEMENTS[key] * value
            for key, value in self.composition.items()
        )

    @property
    def phase(self):
        return PHASE_NUM[self._phase]

    @phase.setter
    def phase(self, value):
        if value is None:
            return
        self._phase = value.strip()

    @property
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self, value):
        self._symbol = value.strip()


class ThermoContainer(BaseDictContainer):
    """Thermo dict container object (storing datas)."""
    _store_type = Thermo
