################################################################################
# @file mechanism.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
from copy import deepcopy
from .base import Base
from .element import ElementContainer
from .species import SpeciesContainer
from .reaction import ReactionContainer, def_unit_k0, def_unit_Ea
from .thermo import ThermoContainer
from .transport import TransportContainer


class Mechanism(Base):
    """Object for storing mechanism data."""

    def __add__(self, other):
        if type(self) == type(other):
            if other is None:
                return deepcopy(self)

                new = deepcopy(self)
                new._combine(other)
                return new
        else:
            raise(TypeError(type(self), type(other)))

    def __init__(self, name, elements=None, specii=None,
                 reactions=None, thermos=None, transports=None):
        self.name = name
        self.elements = elements
        self.specii = specii
        self.reactions = reactions
        self.thermos = thermos
        self.transports = transports

    def __str__(self):
        return self.name

    def _check_list(self):
        reaction_strings = [x.reaction_string_short() for x in self.reactions]

        return [
            (self.elements.is_valid(), 'Element data not valid!'),
            (self.specii.is_valid(), 'Species data not valid!'),
            (self.reactions.is_valid(elements=self.elements, specii=self.specii,
             reaction_strings=reaction_strings), 'Reaction data not valid!'),
            (self.thermos.is_valid(), 'Thermo data not valid!'),
            (self.transports.is_valid(), 'Transport data not valid!'),
            (len(self.inert_specii()) >= 1, 'Missing inert species!'),
        ]

    def _combine(self, other):
        for key, value in other.elements.items():
            self.elements[key] = value

        for key, value in other.thermos.items():
            self.thermos[key] = value

        for key, value in other.transports.items():
            self.transports[key] = value

        for key, value in other.specii.items():
            self.specii[key] = value

        for value in other.reactions:
            self.reactions.append(value)

        self._update_thermos()
        self._update_transports()

    def _update_thermos(self):
        for key in self.thermos:
            if key in self.specii:
                self.specii[key].thermo = self.thermos[key]

    def _update_transports(self):
        for key in self.transports:
            if key in self.specii:
                self.specii[key].transport = self.transports[key]

    def chemkinify(self, prefix, unit_k0=None, unit_Ea=None):
        with open('{:}.mech'.format(prefix), 'w') as fp:
            fp.writelines(
                ['ELEMENTS\n'] + self.elements.chemkinify() + ['END\n']
            )
            fp.writelines(
                ['SPECIES\n'] + self.specii.chemkinify() + ['END\n']
            )

            unit_k0 = self.reactions.unit_k0 if unit_k0 is None else unit_k0
            unit_Ea = self.reactions.unit_Ea if unit_Ea is None else unit_Ea

            tmp = 'REACTIONS'
            if unit_k0 != def_unit_k0:
                tmp += ' {:}'.format(unit_k0)
            if unit_Ea != def_unit_Ea:
                tmp += ' {:}'.format(unit_Ea)

            fp.writelines(
                ['{:}\n'.format(tmp)] + self.reactions.chemkinify(
                    unit_k0=unit_k0, unit_Ea=unit_Ea) + ['END\n']
            )

        with open('{:}.thermo'.format(prefix), 'w') as fp:
            fp.writelines([
                'THERMO\n',
                '{:}\n'.format(' '.join(str(x) for x in self.thermos.bounds))] +
                self.thermos.chemkinify(**{'keys': self.specii.keys()}) +
                ['END\n']
            )

        with open('{:}.transport'.format(prefix), 'w') as fp:
            fp.writelines(self.transports.chemkinify(
                **{'keys': self.specii.keys()}))

    def inert_specii(self):
        species_list = [y for x in self.reactions for y in x.reactants]
        species_list += [y for x in self.reactions for y in x.products]
        return [x for x in self.specii.keys()
                if x not in list(set(species_list))]

    def overwrite_all_flags(self, state):
        self.elements.all_flag = state
        self.specii.all_flag = state
        self.reactions.all_flag = state
        self.thermos.all_flag = state
        self.transports.all_flag = state

    @property
    def thermos(self):
        return self._thermos

    @thermos.setter
    def thermos(self, value):
        self._thermos = value
        self._update_thermos()

    @property
    def transports(self):
        return self._transports

    @transports.setter
    def transports(self, value):
        self._transports = value
        self._update_transports()
