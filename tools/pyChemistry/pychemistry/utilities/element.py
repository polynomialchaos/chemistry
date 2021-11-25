################################################################################
# @file element.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
from .base import Base, BaseDictContainer
from .constants import REF_ELEMENTS
from .utilities import chunk_list


class Element(Base):
    """Element object (storing data)."""

    def __init__(self, symbol, mass=None):
        self.symbol = symbol
        self.mass = mass

    def _data_check(self):
        return [
            (self.mass > 0.0,
             'Atomic mass not valid (am={:})'.format(self.mass)),
        ]

    def chemkinify(self):
        tmp = '/{:}/'.format(self._mass) if hasattr(self, '_mass') else ''
        return self.symbol + tmp

    @property
    def mass(self):
        return REF_ELEMENTS[self.symbol]

    @mass.setter
    def mass(self, value):
        if value is None:
            return
        self._mass = float(value)
        REF_ELEMENTS[self.symbol] = self._mass

    @property
    def symbol(self):
        return self._symbol.upper()

    @symbol.setter
    def symbol(self, value):
        self._symbol = value.strip()


class ElementContainer(BaseDictContainer):
    """Element dict container object (storing datas)."""
    _store_type = Element

    def __str__(self):
        return self.symbol

    def chemkinify(self, keys=None, lineLength=80, **kwargs):
        tmp = self.keys() if keys is None else keys
        element_strings = [self.__getitem__(
            key).chemkinify(**kwargs) for key in tmp]
        max_len = max(len(x) for x in element_strings) + 1

        return [
            ' '.join(['{1:{0:}}'.format(max_len, x) for x in chunk]) + '\n'
            for chunk in chunk_list(element_strings, int(lineLength / max_len))
        ]
