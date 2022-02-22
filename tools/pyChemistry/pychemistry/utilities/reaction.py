################################################################################
# @file reaction.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import logging
from enum import Enum, unique
from .base import Base, BaseListContainer
from .utilities import as_short, chunk_list
from .constants import NA, CAL_JOULE, RM, CM_M


@unique
class ReactionType(Enum):
    """Reaction type enumeration."""
    DEFAULT = 0  # Reaction type identifier for default reactions
    THREE_BODY = 1  # Reaction type identifier for three-body reactions
    PRESSURE = 2  # Reaction type identifier for pressure dependent reactions


conv_k0 = {'MOLES': 1.0, 'MOLECULES': 1 / NA}
conv_Ea = {'CAL/MOLE': CAL_JOULE, 'KCAL/MOLE': CAL_JOULE * 1000.0,
           'JOULES/MOLE': 1.0, 'KJOULES/MOLE': 1000.0, 'KELVINS': RM}
DEF_UNIT_K0 = 'MOLES'
DEF_UNIT_EA = 'CAL/MOLE'


class Reaction(Base):
    """Reaction object (storing data)."""
    # pylint: disable=too-many-instance-attributes
    # pylint: disable=too-many-public-methods

    def __init__(self, reaction_type, reactants, products, is_reversible):
        self.reaction_type = reaction_type
        self.reactants = reactants
        self.products = products
        self.is_reversible = is_reversible

    def _add_species_string(self):
        if self.reaction_type == ReactionType.THREE_BODY:
            return '+M'

        if self.reaction_type == ReactionType.PRESSURE:
            return '(+{:})'.format(self.falloff_species)

        return ''

    def _check_list(self, **kwargs):
        elements = kwargs['elements']
        specii = kwargs['specii']
        reaction_strings = kwargs['reaction_strings']

        string = self.reaction_string_short()
        n_duplicates = reaction_strings.count(string)
        is_duplicate = self.is_duplicate()
        last_troe = (self.troe_coeff[3] != 0.0) if self.troe_coeff and len(
            self.troe_coeff) == 4 else True

        return [
            (len(self.reactants.keys()) > 0, 'No reactants are speciefied!'),
            (len(self.products.keys()) > 0, 'No products are speciefied!'),
            (last_troe,
             'Last TROE coefficient must not be zero (set to SMALL)!'),
            (not (n_duplicates == 1 and is_duplicate),
             'Reaction has duplicate flag, but no duplicate reaction is found!'),
            (not (n_duplicates > 1 and not is_duplicate),
             'Reaction is duplicate, but no duplicate flag is provided!'),
            (check_element_balance(elements, specii, self),
             'Reaction does not balance in elements!'),
            (check_mass_balance(specii, self),
             'Reaction does not balance in mass!'),
        ]

    def __str__(self):
        return self.reaction_string()

    def _species_string(self, values):
        tmp = '+'.join(values)
        return tmp + self._add_species_string()

    def _stoich_species_string(self, values):
        tmp = '+'.join(
            ['{:}{:}'.format(as_short(values[key],
                                      round_digits=12), key) for key in values]
        )
        return tmp + self._add_species_string()

    @property
    def adv_arr_key(self):
        """Reaction (HIGH/LOW) key."""
        return self._adv_arr_key

    @adv_arr_key.setter
    def adv_arr_key(self, value):
        self._adv_arr_key = value

    @property
    def adv_arr_coeff(self):
        """Reaction (HIGH/LOW) Arrhenius coefficients."""
        return self._adv_arr_coeff if self.has_adv_arr_coeff else []

    @adv_arr_coeff.setter
    def adv_arr_coeff(self, value):
        self._adv_arr_coeff = [float(x) for x in value]

    @property
    def arr_coeff(self):
        """Reaction Arrhenius coefficients."""
        return self._arr_coeff

    @arr_coeff.setter
    def arr_coeff(self, value):
        self._arr_coeff = [float(x) for x in value]

    def chemkinify(self, **kwargs):
        min_len = kwargs.get('min_len', 45)
        reaction_string = kwargs.get('reaction_string', self.reaction_string())
        unit_k0 = kwargs.get('unit_k0', DEF_UNIT_K0)
        unit_ea = kwargs.get('unit_ea', DEF_UNIT_EA)

        reaction_strings = []
        reaction_strings.append(
            '{1:{0:}} {2:10.3e} {3:10.4f} {4:10.3f}'.format(
                min_len, reaction_string,
                *si_to_cmsk(self.arr_coeff, self.f_conv_si(), unit_k0, unit_ea))
        )

        if self.has_efficiencies:
            efficiencies = [
                '{:}/{:.2f}/'.format(key, self.efficiencies[key])
                for key in self.efficiencies
            ]
            reaction_strings.extend(
                ['  ' + ' '.join(x) for x in chunk_list(efficiencies, 10)]
            )

        if self.has_rev_arr_coeff:
            reaction_strings.append(
                '  REV /{:10.3e} {:10.4f} {:10.3f}/'.format(*si_to_cmsk(
                    self.rev_arr_coeff, self.r_conv_si(), unit_k0, unit_ea))
            )

        if self.has_adv_arr_coeff and self.adv_arr_key == 'LOW':
            reaction_strings.append(
                '  LOW /{:10.3e} {:10.4f} {:10.3f}/'.format(*si_to_cmsk(
                    self.adv_arr_coeff, self.f_conv_si(add=1.0),
                    unit_k0, unit_ea))
            )
        elif self.has_adv_arr_coeff and self.adv_arr_key == 'HIGH':
            reaction_strings.append(
                '  HIGH/{:10.3e} {:10.4f} {:10.3f}/'.format(*si_to_cmsk(
                    self.adv_arr_coeff, self.f_conv_si(add=-1.0),
                    unit_k0, unit_ea))
            )

        if self.has_troe_coeff:
            reaction_strings.append(
                '  TROE/{:}/'.format(' '.join(['{:.4e}'.format(x)
                                               for x in self.troe_coeff]))
            )

        if self.has_f_orders:
            reaction_strings.extend(
                '  FORD/{:} {:.6f}/'.format(key, self.f_orders[key])
                for key in self.f_orders
            )

        if self.has_r_orders:
            reaction_strings.extend(
                '  RORD/{:} {:.6f}/'.format(key, self.r_orders[key])
                for key in self.r_orders
            )

        if self.has_flags:
            reaction_strings.extend(
                '  {:}'.format(key) for key in self.flags
            )

        return '\n'.join(reaction_strings)

    @property
    def efficiencies(self):
        """Reaction third-body efficiencies."""
        return self._efficiencies if self.has_efficiencies else {}

    @efficiencies.setter
    def efficiencies(self, value):
        self._efficiencies = {key: float(value)
                              for key, value in value.items()}

    @property
    def falloff_species(self):
        """Reaction fall-off species (three-body or pressure)."""
        return self._falloff_species

    @falloff_species.setter
    def falloff_species(self, value):
        self._falloff_species = value.strip()

    @property
    def flags(self):
        """Reaction auxiliary flags."""
        return self._flags if self.has_flags else []

    @flags.setter
    def flags(self, value):
        self._flags = [str(x) for x in value]

    def f_conv_si(self, add=0.0):
        """Reaction forward conversion factor."""
        sum_nu = self.f_order() + add
        return CM_M ** (3 * (sum_nu - 1))

    def f_order(self):
        """Reaction forward order."""
        for_order = sum(self.reactant_orders)
        for_order += (1.0 if self.reaction_type ==
                      ReactionType.THREE_BODY else 0.0)
        return for_order

    @property
    def f_orders(self):
        """Reaction forward orders."""
        return self._f_orders if self.has_f_orders else {}

    @f_orders.setter
    def f_orders(self, value):
        self._f_orders = {key: float(value) for key, value in value.items()}

    def is_duplicate(self):
        """Reaction duplicate flag."""
        return any('DUPL' in x for x in self.flags)

    @property
    def has_adv_arr_coeff(self):
        """Reaction has (HIGH/LOW) Arrhenius coefficients."""
        return hasattr(self, '_adv_arr_coeff')

    @property
    def has_efficiencies(self):
        """Reaction has third-body efficiencies."""
        return hasattr(self, '_efficiencies')

    @property
    def has_flags(self):
        """Reaction has auxiliary flags."""
        return hasattr(self, '_flags')

    @property
    def has_f_orders(self):
        """Reaction has forward orders."""
        return hasattr(self, '_f_orders')

    @property
    def has_r_orders(self):
        """Reaction has backward (reverse) orders."""
        return hasattr(self, '_r_orders')

    @property
    def has_rev_arr_coeff(self):
        """Reaction has reverse Arrhenius coefficients."""
        return hasattr(self, '_rev_arr_coeff')

    @property
    def has_troe_coeff(self):
        """Reaction has TROE coefficients."""
        return hasattr(self, '_troe_coeff')

    @property
    def is_reversible(self):
        """Reaction reversible flag."""
        return self._is_reversible

    @is_reversible.setter
    def is_reversible(self, value):
        self._is_reversible = value

    @property
    def product_orders(self):
        """Reaction product orders."""
        return [self.r_orders.get(key, self.products[key])
                for key in self.products]

    @property
    def products(self):
        """Reaction products."""
        return self._products

    @products.setter
    def products(self, value):
        self._products = {key: float(value) for key, value in value.items()}

    def r_conv_si(self, add=0.0):
        """Reaction backward (reverse) conversion factor."""
        sum_nu = self.r_order() + add
        return CM_M ** (3 * (sum_nu - 1))

    def r_order(self):
        """Reaction backward (reverse) order."""
        rev_order = sum(self.product_orders)
        rev_order += (1.0 if self.reaction_type ==
                      ReactionType.THREE_BODY else 0.0)
        return rev_order

    @property
    def r_orders(self):
        """Reaction backward (reverse) orders."""
        return self._r_orders if self.has_r_orders else {}

    @r_orders.setter
    def r_orders(self, value):
        self._r_orders = {key: float(value) for key, value in value.items()}

    @property
    def reactant_orders(self):
        """Reaction reactant orders."""
        return [self.f_orders.get(key, self.reactants[key])
                for key in self.reactants]

    @property
    def reactants(self):
        """Reaction reactants."""
        return self._reactants

    @reactants.setter
    def reactants(self, value):
        self._reactants = {key: float(value) for key, value in value.items()}

    def reaction_string(self):
        """Reaction string (including stoich. values)."""
        delimiter = '<=>' if self.is_reversible else '=>'
        return ''.join([self._stoich_species_string(self.reactants), delimiter,
                        self._stoich_species_string(self.products)])

    def reaction_string_short(self):
        """Reaction short string (excluding stoich. values)."""
        delimiter = '<=>' if self.is_reversible else '=>'
        return ''.join([self._species_string(self.reactants), delimiter,
                        self._species_string(self.products)])

    @property
    def reaction_type(self):
        """Reaction type."""
        return self._reaction_type

    @reaction_type.setter
    def reaction_type(self, value):
        self._reaction_type = value

    @property
    def rev_arr_coeff(self):
        """Reaction reverse Arrhenius coefficients."""
        return self._rev_arr_coeff if self.has_rev_arr_coeff else []

    @rev_arr_coeff.setter
    def rev_arr_coeff(self, value):
        self._rev_arr_coeff = [float(x) for x in value]

    @property
    def sum_nu(self):
        """Sum of stoichiometric values (products - educts)."""
        return sum(self.products.values()) - sum(self.reactants.values())

    @property
    def troe_coeff(self):
        """Reaction TROE coefficients."""
        return self._troe_coeff if self.has_troe_coeff else []

    @troe_coeff.setter
    def troe_coeff(self, value):
        self._troe_coeff = [float(x) for x in value]


class ReactionContainer(BaseListContainer):
    """Reaction list container object (storing datas)."""
    _store_type = Reaction

    def __init__(self, items=None, all_flag=False):
        super().__init__(items=items, all_flag=all_flag)
        self.unit_k0 = DEF_UNIT_K0
        self.unit_ea = DEF_UNIT_EA

    def chemkinify(self, **kwargs):
        unit_k0 = kwargs.get('unit_k0', DEF_UNIT_K0)
        unit_ea = kwargs.get('unit_ea', DEF_UNIT_EA)

        reaction_strings = [x.reaction_string() for x in self]
        max_len = max(len(x) for x in reaction_strings) + 1

        return [
            '! --- Reaction #{:} ---\n{:}\n'.format(
                idx + 1, x.chemkinify(
                    min_len=max_len, reaction_string=reaction_strings[idx],
                    unit_k0=unit_k0, unit_ea=unit_ea)
            ) for idx, x in enumerate(self)
        ]


def cmsk_to_si(coeff, si_conv, unit_k0, unit_ea):
    """Convert the given coefficients to the specified units."""
    return [
        coeff[0] * si_conv * conv_k0[unit_k0],
        coeff[1],
        coeff[2] * conv_Ea[unit_ea],
    ]


def si_to_cmsk(coeff, si_conv, unit_k0, unit_ea):
    """Convert the given coefficients to the specified units."""
    return [
        coeff[0] / si_conv / conv_k0[unit_k0],
        coeff[1],
        coeff[2] / conv_Ea[unit_ea],
    ]


def check_element_balance(elements, species, reaction, limit=1e-5):
    """Check the given reaction for their element balance."""
    elm_sum = {key: 0.0 for key in elements.keys()}

    for spc, spc_nu in reaction.reactants.items():
        for key in species[spc].thermo.composition:
            elm_sum[key] += species[spc].thermo.composition[key] * spc_nu

    for spc, spc_nu in reaction.products.items():
        for key in species[spc].thermo.composition:
            elm_sum[key] -= species[spc].thermo.composition[key] * spc_nu

    if any(abs(x) > limit for x in elm_sum.values()):
        logging.warning('Element conservation violation: %s!', elm_sum)

    return all(abs(x) <= limit for x in elm_sum.values())


def check_mass_balance(species, reaction, limit=1e-5):
    """Check the given reaction for their mass balance."""
    mass_sum = 0.0

    for spc, spc_nu in reaction.reactants.items():
        mass_sum += species[spc].molar_mass * spc_nu

    for spc, spc_nu in reaction.products.items():
        mass_sum -= species[spc].molar_mass * spc_nu

    if abs(mass_sum) > limit:
        logging.warning('Mass conservation violation: %s!', mass_sum)

    return abs(mass_sum) <= limit
