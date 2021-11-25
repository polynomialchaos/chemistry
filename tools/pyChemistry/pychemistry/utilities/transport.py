################################################################################
# @file transport.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
from .base import Base, BaseDictContainer
from .constants import ANGSTROM_SI, DEBYE_SI


CGS_TO_SI_COLLJ = ANGSTROM_SI   # Conversion: Lennard-Jones collision diametre
CGS_TO_SI_DIPMO = DEBYE_SI      # Conversion: Dipole monent
CGS_TO_SI_POL = ANGSTROM_SI**3  # Conversion: Polarizeability


class Transport(Base):
    """Transport object (storing data)."""

    def __init__(self, symbol, geom=None, pot_lj=None, col_lj=None,
                 dip_mo=None, pol=None, rot_rel=None):
        self.symbol = symbol
        self.geom = geom
        self.pot_lj = pot_lj
        self.col_lj = col_lj
        self.dip_mo = dip_mo
        self.pol = pol
        self.rot_rel = rot_rel

    def __str__(self):
        return self.symbol

    def _check_list(self):
        return [
            (self.geom in [0, 1, 2], 'Geometrical configuration not valid!'),
            (self.pot_lj >= 0.0, 'Lennard-Jones potential not valid!'),
            (self.col_lj >= 0.0, 'Lennard-Jones collision diameter not valid!'),
            (self.dip_mo >= 0.0, 'Dipole moment not valid!'),
            (self.pol >= 0.0, 'Polarizability not valid!'),
            (self.rot_rel >= 0.0, 'Rotational relaxation not valid!'),
        ]

    def chemkinify(self):
        return '{:<16}{:<4}{:10.3f} {:10.3f} {:10.3f} {:10.3f}{:10.3f}'.format(
            self.symbol, self.geom, self.pot_lj,
            self.col_lj / CGS_TO_SI_COLLJ,
            self.dip_mo / CGS_TO_SI_DIPMO,
            self.pol / CGS_TO_SI_POL,
            self.rot_rel
        )

    @property
    def col_lj(self):
        return self._col_lj

    @col_lj.setter
    def col_lj(self, value):
        if value is None:
            return
        self._col_lj = float(value)

    @property
    def dip_mo(self):
        return self._dip_mo

    @dip_mo.setter
    def dip_mo(self, value):
        if value is None:
            return
        self._dip_mo = float(value)

    @property
    def geom(self):
        return self._geom

    @geom.setter
    def geom(self, value):
        if value is None:
            return
        self._geom = int(value)

    @property
    def pol(self):
        return self._pol

    @pol.setter
    def pol(self, value):
        if value is None:
            return
        self._pol = float(value)

    @property
    def pot_lj(self):
        return self._pot_lj

    @pot_lj.setter
    def pot_lj(self, value):
        if value is None:
            return
        self._pot_lj = float(value)

    @property
    def rot_rel(self):
        return self._rot_rel

    @rot_rel.setter
    def rot_rel(self, value):
        if value is None:
            return
        self._rot_rel = float(value)

    @property
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self, value):
        self._symbol = value.strip()


class TransportContainer(BaseDictContainer):
    """Transport dict container object (storing datas)."""
    _store_type = Transport
