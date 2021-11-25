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


class Species(Base):
    """Species object (storing data)."""

    def __init__(self, symbol, molar_mass=None, thermo=None, transport=None):
        self.symbol = symbol
        self.molar_mass = molar_mass
        self.thermo = thermo
        self.transport = transport

    def __str__(self):
        return self.symbol

    def _check_list(self):
        return [
            (self.molar_mass >= 0.0,
             'Molar mass not valid (mm={:})!'.format(self.molar_mass)),
            (self.thermo is not None, 'Missing thermo data!'),
            (self.transport is not None, 'Missing transport data!'),
        ]

    def chemkinify(self):
        return self.symbol

    @property
    def molar_mass(self):
        if hasattr(self, '_molar_mass'):
            return self._molar_mass
        elif hasattr(self, '_thermo') and self.thermo:
            return self.thermo.molar_mass
        else:
            raise(ValueError(
                'Neither molar_mass nor thermo data for {:} specified!'.format(
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
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self, value):
        self._symbol = value.strip()

    @property
    def thermo(self):
        return self._thermo

    @thermo.setter
    def thermo(self, value):
        if value is None:
            return
        self._thermo = value

    @property
    def transport(self):
        return self._transport

    @transport.setter
    def transport(self, value):
        if value is None:
            return
        self._transport = value


class SpeciesContainer(ElementContainer):
    """Species dict container object (storing datas)."""
    _store_type = Species
