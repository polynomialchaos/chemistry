####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import numpy as np
from collections import OrderedDict

from .base import Base, BaseOrderedDictContainer
from .constants import REF_ELEMENTS, RM

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
PHASE_NUM = {'G': 0}  # Conversion: Phase string
N_BOUNDS = 3         # Number of temperature bounds
N_NASA = 7         # Number of NASA polynomial coefficients

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
class Thermo(Base):
    """Object for storing thermo data. Inputs in SI units."""

    def __init__(self, symbol, info='', composition=None, phase='',
                 bounds=None, coeff_low=None, coeff_high=None):

        self.symbol = symbol

        self.info = info
        self.composition = OrderedDict() if composition is None else composition
        self.phase = phase
        self.bounds = [] if bounds is None else bounds
        self.coeff_low = [] if coeff_low is None else coeff_low
        self.coeff_high = [] if coeff_high is None else coeff_high

    def __str__(self):
        return self.symbol

    def _checklist(self):
        return [
            (bool(self.composition), 'Elemental composition empty'),
            (len(self.bounds) == N_BOUNDS,
             'NASA temp. bounds of wrong size ({:})'.format(len(self.bounds))),
            (len(self.coeff_low) == N_NASA,
             'NASA low temp. coefficients of wrong size ({:})'.format(len(self.coeff_low))),
            (len(self.coeff_high) == N_NASA,
             'NASA high temp. coefficients of wrong size ({:})'.format(len(self.coeff_high))),
            (self.phase == 0, 'Phase not valid'),
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
        tmp = sorted([float(x) for x in value])
        self._bounds = [min(tmp), max(tmp)]
        self._bounds += [x for x in tmp if x not in self._bounds]

    def chemkinify(self):
        symb_info = '{:} {:}'.format(self.symbol, self.info)
        short_comps, long_comps = self._get_comp_strings()
        return '{:79}1{:}\n{:79}2\n{:79}3\n{:79}4'.format(
            '{:24}{:20}{:1}{:-10.3f}{:-10.3f}{:-10.3f}'.format(
                symb_info[:24], short_comps, self._phase, *self._bounds
            ),
            '' if not long_comps else ' &\n{:}'.format(long_comps),
            ''.join(['{: 15.8E}'.format(coeff)
                    for coeff in self.coeff_high[:5]]),
            ''.join(['{: 15.8E}'.format(coeff)
                    for coeff in self.coeff_high[5:] + self.coeff_low[:3]]),
            ''.join(['{: 15.8E}'.format(coeff)
                    for coeff in self.coeff_low[3:]]),
        )

    @property
    def coeff_high(self):
        return self._coeff_high

    @coeff_high.setter
    def coeff_high(self, value):
        self._coeff_high = [float(x) for x in value]

    @property
    def coeff_low(self):
        return self._coeff_low

    @coeff_low.setter
    def coeff_low(self, value):
        self._coeff_low = [float(x) for x in value]

    @property
    def composition(self):
        return self._composition

    @composition.setter
    def composition(self, value):
        self._composition = OrderedDict(
            (key, float(value[key])) for key in value
        )

    def fit_nasa_polynomial(self, ref_T, cp, h298, s298):
        polyFit = np.flip(np.polyfit(ref_T, cp, deg=(7-1-2)), axis=0)
        tmpCoeff = [x for x in polyFit] + [0.0, 0.0]

        tmpCoeff[-2:] = [
            h298 - calc_dimless_h(298.0, tmpCoeff) * 298.0,
            s298 - calc_dimless_s(298.0, tmpCoeff)
        ]

        self.coeff_low = [x * self.molar_mass / RM for x in tmpCoeff]
        self.coeff_high = [x for x in self.coeff_low]

    @property
    def info(self):
        return self._info

    @info.setter
    def info(self, value):
        self._info = value.strip()

    @property
    def molar_mass(self):
        return sum(
            REF_ELEMENTS[key] * self.composition[key]
            for key in self.composition
        )

    @property
    def phase(self):
        return PHASE_NUM[self._phase]

    @phase.setter
    def phase(self, value):
        self._phase = value.strip()

    @property
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self, value):
        self._symbol = value.strip()

class ThermoContainer(BaseOrderedDictContainer):
    """Container (storage) object for Thermo class objects."""
    _type = Thermo

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
def calc_dimless_cp(T, a):
    """Return the dimensionless heat capacity in constant pressure polynomial at the given temperature(s)."""
    return a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4])))

def calc_dimless_h(T, a):
    """Return the dimensionless enthalpy polynomial at the given temperature(s)."""
    return a[0] + T * (a[1] / 2 + T * (a[2] / 3 + T * (a[3] / 4 + T * a[4] / 5))) + a[5] / T

def calc_dimless_s(T, a):
    """Return the dimensionless entropy polynomial at the given temperature(s)."""
    return a[0] * np.log(T) + T * (a[1] + T * (a[2] / 2 + T * (a[3] / 3 + T * a[4] / 4))) + a[6]

def calc_dimless_g(T, a):
    """Return the dimensionless free gibbs energy polynomial at the given temperature(s)."""
    return a[0] * (1 - np.log(T)) - T * (a[1] / 2 + T * (a[2] / 6 + T * (a[3] / 12 + T * a[4] / 20))) + a[5] / T - a[6]
