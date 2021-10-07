####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
from .base import Base, BaseOrderedDictContainer
from .constants import ANGSTROM_SI, DEBYE_SI

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
CGS_TO_SI_COLLJ = ANGSTROM_SI       # Conversion: Lennard-Jones collision diametre
CGS_TO_SI_DIPMO = DEBYE_SI          # Conversion: Dipole monent
CGS_TO_SI_POL = ANGSTROM_SI**3    # Conversion: Polarizeability

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
class Transport(Base):
    """Object for storing transport data. Inputs in SI units."""

    def __init__(self, symbol, geom=0, pot_lj=0.0, col_lj=0.0, dip_mo=0.0, pol=0.0, rot_rel=0.0):

        self.symbol = symbol

        self.geom = geom
        self.pot_lj = pot_lj
        self.col_lj = col_lj
        self.dip_mo = dip_mo
        self.pol = pol
        self.rot_rel = rot_rel

    def __str__(self):
        return self.symbol

    def _checklist(self):
        return [
            (self.geom in [0, 1, 2], 'Geometrical configuration not valid'),
            (self.pot_lj >= 0.0, 'Lennard-Jones potential not valid'),
            (self.col_lj >= 0.0, 'Lennard-Jones collision diameter not valid!'),
            (self.dip_mo >= 0.0, 'Dipole moment not valid!'),
            (self.pol >= 0.0, 'Polarizability not valid!'),
            (self.rot_rel >= 0.0, 'Rotational relaxation not valid'),
        ]

    def chemkinify(self):
        return '{:<16}{:<4}{:10.3f} {:10.3f} {:10.3f} {:10.3f}{:10.3f}'.format(
            self.symbol, self.geom, self.pot_lj, self.col_lj / CGS_TO_SI_COLLJ,
            self.dip_mo / CGS_TO_SI_DIPMO, self.pol / CGS_TO_SI_POL, self.rot_rel
        )

    @property
    def col_lj(self):
        return self._col_lj

    @col_lj.setter
    def col_lj(self, value):
        self._col_lj = float(value)

    @property
    def dip_mo(self):
        return self._dip_mo

    @dip_mo.setter
    def dip_mo(self, value):
        self._dip_mo = float(value)

    @property
    def geom(self):
        return self._geom

    @geom.setter
    def geom(self, value):
        self._geom = int(value)

    @property
    def pol(self):
        return self._pol

    @pol.setter
    def pol(self, value):
        self._pol = float(value)

    @property
    def pot_lj(self):
        return self._pot_lj

    @pot_lj.setter
    def pot_lj(self, value):
        self._pot_lj = float(value)

    @property
    def rot_rel(self):
        return self._rot_rel

    @rot_rel.setter
    def rot_rel(self, value):
        self._rot_rel = float(value)

    @property
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self, value):
        self._symbol = value.strip()

class TransportContainer(BaseOrderedDictContainer):
    """Container (storage) object for Transport class objects."""
    _type = Transport

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
