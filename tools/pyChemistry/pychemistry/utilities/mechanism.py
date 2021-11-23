################################################################################
# @file mechanism.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import sys
from collections import OrderedDict
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
        if other is None:
            return deepcopy(self)

        if isinstance(self, other.__class__):
            new = deepcopy(self)

            new.name += '+{:}'.format(other.name)
            new.models.update({key: other.models[key] for key in other.models})

            for key in other.elements.keys():
                new.elements[key] = other.elements[key]

            for key in other.specii.keys():
                new.specii[key] = other.specii[key]

            for x in other.reactions:
                new.reactions.append(x)

            for key in other.thermos.keys():
                new.thermos[key] = other.thermos[key]

            for key in other.transports.keys():
                new.transports[key] = other.transports[key]

            new._update_thermos()
            new._update_transports()

            return new
        else:
            raise(TypeError(Mechanism, type(other)))

    def __init__(self, name, elements=None, specii=None,
                 reactions=None, thermos=None, transports=None, models=None):

        self.name = name
        self.elements = ElementContainer() if elements is None else elements
        self.specii = SpeciesContainer() if specii is None else specii
        self.reactions = ReactionContainer() if reactions is None else reactions
        self.thermos = ThermoContainer() if thermos is None else thermos
        self.transports = TransportContainer() \
            if transports is None else transports
        self.models = OrderedDict() if models is None else models

    def __str__(self):
        return self.name

    def _checklist(self):
        reaction_strings = [x.reaction_string_short() for x in self.reactions]

        return [
            (self.elements.is_valid(), 'Element data not valid'),
            (self.specii.is_valid(), 'Species data not valid'),
            (self.reactions.is_valid(
                elements=self.elements, specii=self.specii,
                reaction_strings=reaction_strings),
                'Reaction data not valid'
             ),
            (self.thermos.is_valid(), 'Thermo data not valid'),
            (self.transports.is_valid(), 'Transport data not valid'),
            (len(self.inert_specii()) >= 1, 'Missing inert species'),
        ]

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
            fp.writelines(['ELEMENTS\n'] +
                          self.elements.chemkinify() + ['END\n'])
            fp.writelines(['SPECIES\n'] + self.specii.chemkinify() + ['END\n'])

            if unit_k0 is None:
                unit_k0 = self.reactions.unit_k0
            if unit_Ea is None:
                unit_Ea = self.reactions.unit_Ea

            tmp = 'REACTIONS'
            if unit_k0 != def_unit_k0:
                tmp += ' {:}'.format(unit_k0)
            if unit_Ea != def_unit_Ea:
                tmp += ' {:}'.format(unit_Ea)

            fp.writelines(['{:}\n'.format(tmp)] +
                          self.reactions.chemkinify(
                unit_k0=unit_k0, unit_Ea=unit_Ea) + ['END\n'])

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
