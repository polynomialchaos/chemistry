################################################################################
# @file mechanism.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
from copy import deepcopy
from .base import Base
from .reaction import DEF_UNIT_K0, DEF_UNIT_EA


class Mechanism(Base):
    """Mechanism object (storing data)."""
    # pylint: disable=too-many-instance-attributes

    def __add__(self, other):
        if isinstance(other, type(self)):
            if other is None:
                return deepcopy(self)

            new = deepcopy(self)
            new._combine(other)
            return new

        raise TypeError('Wrong type "{:}" ({:})!'.format(
            type(other), type(self)))

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

    def _check_list(self, **_):
        reaction_strings = [x.reaction_string_short() for x in self.reactions]

        return [
            (self.elements.is_valid(), 'Element data not valid!'),
            (self.specii.is_valid(), 'Species data not valid!'),
            (self.reactions.is_valid(elements=self.elements, specii=self.specii,
                                     reaction_strings=reaction_strings),
             'Reaction data not valid!'),
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

    def chemkinify(self, **kwargs):
        """Return a CHEMKIN formatted string."""
        prefix = kwargs['prefix']
        unit_k0 = kwargs.get('unit_k0', self.reactions.unit_k0)
        unit_ea = kwargs.get('unit_ea', self.reactions.unit_ea)

        with open('{:}.mech'.format(prefix), 'w') as fptr:
            fptr.writelines(
                ['ELEMENTS\n'] + self.elements.chemkinify() + ['END\n']
            )
            fptr.writelines(
                ['SPECIES\n'] + self.specii.chemkinify() + ['END\n']
            )

            tmp = 'REACTIONS'
            if unit_k0 != DEF_UNIT_K0:
                tmp += ' {:}'.format(unit_k0)
            if unit_ea != DEF_UNIT_EA:
                tmp += ' {:}'.format(unit_ea)

            fptr.writelines(
                ['{:}\n'.format(tmp)] + self.reactions.chemkinify(
                    unit_k0=unit_k0, unit_ea=unit_ea) + ['END\n']
            )

        with open('{:}.thermo'.format(prefix), 'w') as fptr:
            fptr.writelines([
                'THERMO\n',
                '{:}\n'.format(' '.join(str(x) for x in self.thermos.bounds))] +
                self.thermos.chemkinify(**{'keys': self.specii.keys()}) +
                ['END\n']
            )

        with open('{:}.transport'.format(prefix), 'w') as fptr:
            fptr.writelines(self.transports.chemkinify(
                **{'keys': self.specii.keys()}))

    def inert_specii(self):
        """Return a list of inert specii."""
        species_list = [y for x in self.reactions for y in x.reactants]
        species_list += [y for x in self.reactions for y in x.products]
        return [x for x in self.specii.keys()
                if x not in list(set(species_list))]

    @property
    def thermos(self):
        """Thermo datas."""
        return self._thermos

    @thermos.setter
    def thermos(self, value):
        self._thermos = value
        self._update_thermos()

    @property
    def transports(self):
        """Transport datas."""
        return self._transports

    @transports.setter
    def transports(self, value):
        self._transports = value
        self._update_transports()
