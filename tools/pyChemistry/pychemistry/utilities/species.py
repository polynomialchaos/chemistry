################################################################################
# @file species.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import numpy as np
from .base import Base
from .element import ElementContainer
from .constants import NA
from .utilities import chunk_list
from .thermo import calc_dimless_cp, calc_dimless_h
from .thermo import calc_dimless_s, calc_dimless_g


class Species(Base):
    """Object for storing species data. Inputs in SI units."""

    def __init__(self, symbol, molar_mass=None, thermo=None, transport=None):

        self.symbol = symbol

        self.molar_mass = molar_mass
        self.thermo = thermo
        self.transport = transport

    def __str__(self):
        return self.symbol

    def _checklist(self):
        add_checks = []
        if self.thermo:
            T_min, T_max = self.thermo.bounds[0], self.thermo.bounds[2]
            T = np.linspace(T_min, T_max, 100)
            cp_poly = self._get_dimless_cp(T)

            add_checks.append(
                (all(cp_poly >= 0.0),
                    'Cp polynomial not valid (reduce bounds to {:.2f})'.format(
                    T[np.where(cp_poly >= 0)[-1][-1]])),
            )

        return [
            (self.thermo is not None, 'Missing thermo data'),
            (self.transport is not None, 'Missing transport data'),
            (self.molar_mass >= 0.0,
             'Molar mass not valid (mm={:})'.format(self.molar_mass)),
        ] + add_checks

    def _get_dimless_cp(self, T):
        return self._vector_handler(T, calc_dimless_cp)

    def _get_dimless_g(self, T):
        return self._vector_handler(T, calc_dimless_g)

    def _get_dimless_h(self, T):
        return self._vector_handler(T, calc_dimless_h)

    def _get_dimless_s(self, T):
        return self._vector_handler(T, calc_dimless_s)

    def _vector_handler(self, T, function):
        if isinstance(T, (list, np.ndarray)):
            T_l, T_h = splitAt(T, self.thermo.bounds[1])
            return np.concatenate((function(T_l, self.thermo.coeff_low),
                                   function(T_h, self.thermo.coeff_high)))
        else:
            return function(T, self.thermo.coeff_low
                            if T < self.thermo.bounds[1]
                            else self.thermo.coeff_high)

    @property
    def molar_mass(self):
        if hasattr(self, '_molar_mass'):
            return self._molar_mass
        elif self.thermo:
            return self.thermo.molar_mass
        else:
            raise(ValueError('Neither molar mass nor elements/thermo properties '
                             'have been set for species "{:}"!'.format(
                                 self.symbol)))

    @molar_mass.setter
    def molar_mass(self, value):
        if value is None:
            return
        self._molar_mass = float(value)

    @property
    def molecule_weight(self):
        return self.molar_mass / NA

    @property
    def thermo(self):
        return self._thermo

    @thermo.setter
    def thermo(self, x):
        self._thermo = x

    @property
    def transport(self):
        return self._transport

    @transport.setter
    def transport(self, x):
        self._transport = x

    @property
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self, value):
        self._symbol = value.strip()

    def chemkinify(self):
        return self.symbol


class SpeciesContainer(ElementContainer):
    """Container (storage) object for Species class objects."""
    _type = Species


def splitAt(x, position):
    """Return two lists of provided values, splid at the given position."""
    idx = np.where(x > position)[0][0]
    return (x[:idx], x[idx:])
