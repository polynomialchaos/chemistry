####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import logging
from enum import Enum, unique
from collections import OrderedDict

from .base import Base, BaseListContainer
from .utilities import as_short, chunk_list
from .constants import NA, CAL_JOULE, RM, CM_M

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
@unique
class ReactionType( Enum ):
    DEFAULT     = 0 # Reaction type identifier for default reactions
    THREE_BODY  = 1 # Reaction type identifier for three-body reactions
    PRESSURE    = 2 # Reaction type identifier for pressure dependent reactions

conv_k0     = {'MOLES': 1.0, 'MOLECULES': 1 / NA}
conv_Ea     = {'CAL/MOLE': CAL_JOULE, 'KCAL/MOLE': CAL_JOULE * 1000.0, 'JOULES/MOLE': 1.0, 'KJOULES/MOLE': 1000.0, 'KELVINS': RM}
def_unit_k0 = 'MOLES'
def_unit_Ea = 'CAL/MOLE'

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
class Reaction( Base ):
    """Object for storing Reaction data."""
    def __init__( self, reaction_type, reactants, products, is_reversible, arr_coeff=None, falloff_species=None,
        rev_arr_coeff=None, adv_arr_key=None, adv_arr_coeff=None, troe_coeff=None,
        f_orders=None, r_orders=None, flags=None, efficiencies=None, elementar=None ):

        self.reaction_type      = reaction_type
        self.reactants          = reactants
        self.products           = products
        self.is_reversible      = is_reversible

        self.arr_coeff          = [] if arr_coeff is None else arr_coeff
        self.falloff_species    = '' if falloff_species is None else falloff_species

        self.rev_arr_coeff      = [] if rev_arr_coeff is None else rev_arr_coeff
        self.adv_arr_key        = '' if adv_arr_key is None else adv_arr_key
        self.adv_arr_coeff      = [] if adv_arr_coeff is None else adv_arr_coeff
        self.troe_coeff         = [] if troe_coeff is None else troe_coeff
        self.f_orders           = OrderedDict() if f_orders is None else f_orders
        self.r_orders           = OrderedDict() if r_orders is None else r_orders
        self.flags              = [] if flags is None else flags
        self.efficiencies       = OrderedDict() if efficiencies is None else efficiencies

        self.elementar          = elementar

    @property
    def adv_arr_coeff( self ):
        return self._adv_arr_coeff

    @property
    def arr_coeff( self ):
        return self._arr_coeff

    @property
    def efficiencies( self ):
        return self._efficiencies

    @property
    def falloff_species( self ):
        return self._falloff_species

    @property
    def elementar( self ):
        if hasattr( self, '_elementar' ):
            return self._elementar
        else:
            return self.is_elementar()

    def is_duplicate( self ):
        return any( 'DUPL' in x for x in self.flags )

    @property
    def is_reversible( self ):
        return self._is_reversible

    @property
    def product_orders( self ):
        return [self.r_orders.get( key, self.products[key] ) for key in self.products]

    @property
    def products( self ):
        return self._products

    @property
    def reactant_orders( self ):
        return [self.f_orders.get( key, self.reactants[key] ) for key in self.reactants]

    @property
    def reactants( self ):
        return self._reactants

    @property
    def reaction_type( self ):
        return self._reaction_type

    @property
    def rev_arr_coeff( self ):
        return self._rev_arr_coeff

    @property
    def sum_nu( self ):
        return sum( self.products.values() ) - sum( self.reactants.values() )

    @property
    def troe_coeff( self ):
        return self._troe_coeff

    @adv_arr_coeff.setter
    def adv_arr_coeff( self, value ):
        self._adv_arr_coeff = [float( x ) for x in value]

    @arr_coeff.setter
    def arr_coeff( self, value ):
        self._arr_coeff = [float( x ) for x in value]

    @efficiencies.setter
    def efficiencies( self, value ):
        self._efficiencies = OrderedDict(
            (key, float( value[key] )) for key in value
        )

    @falloff_species.setter
    def falloff_species( self, value ):
        self._falloff_species = value.strip()

    @elementar.setter
    def elementar( self, value ):
        if value is None: return
        self._elementar = value

    @is_reversible.setter
    def is_reversible( self, value ):
        self._is_reversible = value

    @products.setter
    def products( self, value ):
        self._products = OrderedDict(
            (key, float( value[key] )) for key in value
        )

    @reactants.setter
    def reactants( self, value ):
        self._reactants = OrderedDict(
            (key, float( value[key] )) for key in value
        )

    @reaction_type.setter
    def reaction_type( self, value ):
        self._reaction_type = value

    @rev_arr_coeff.setter
    def rev_arr_coeff( self, value ):
        self._rev_arr_coeff = [float( x ) for x in value]

    @troe_coeff.setter
    def troe_coeff( self, value ):
        self._troe_coeff = [float( x ) for x in value]

    def is_elementar( self ):
        digit_nu_reactants  = all( self.reactants[key].is_integer() for key in self.reactants )
        digit_nu_products   = all( self.products[key].is_integer() for key in self.products )
        no_ord_key          = not self.f_orders and not self.r_orders
        return bool( digit_nu_reactants and digit_nu_products and no_ord_key )

    def chemkinify( self, min_len=45, reaction_string=None, unit_k0=def_unit_k0, unit_Ea=def_unit_Ea ):
        tmp_string = self.reaction_string() if reaction_string is None else reaction_string
        efficiencies = ['{:}/{:.2f}/'.format( key, self.efficiencies[key] ) for key in self.efficiencies]
        efficiencies = [' '.join( x ) for x in chunk_list( efficiencies, 10 )]

        return '\n'.join(
            (['{1:{0:}} {2:10.3e} {3:10.4f} {4:10.3f}'.format( min_len, tmp_string,
                *si_to_cmsk( self.arr_coeff, self.f_conv_si(), unit_k0, unit_Ea ) )]) +
            (['  {:}'.format( x ) for x in efficiencies]) +
            (['  REV /{:10.3e} {:10.4f} {:10.3f}/'.format(
                *si_to_cmsk( self.rev_arr_coeff, self.r_conv_si(), unit_k0, unit_Ea ) )]
                if self.rev_arr_coeff else []) +
            (['  LOW /{:10.3e} {:10.4f} {:10.3f}/'.format(
                *si_to_cmsk( self.adv_arr_coeff, self.f_conv_si( add=1.0 ), unit_k0, unit_Ea ) )]
                if self.adv_arr_coeff and self.adv_arr_key == 'LOW' else []) +
            (['  HIGH/{:10.3e} {:10.4f} {:10.3f}/'.format(
                *si_to_cmsk( self.adv_arr_coeff, self.f_conv_si( add=-1.0 ), unit_k0, unit_Ea ) )]
                if self.rev_arr_coeff and self.adv_arr_key == 'HIGH' else []) +
            (['  TROE/{:}/'.format( ' '.join( ['{:.4e}'.format( x ) for x in self.troe_coeff] ) )] if self.troe_coeff else []) +
            (['  FORD/{:} {:.6f}/'.format( key, self.f_orders[key] ) for key in self.f_orders]) +
            (['  RORD/{:} {:.6f}/'.format( key, self.r_orders[key] ) for key in self.r_orders]) +
            (['  {:}'.format( key ) for key in self.flags])
        )

    def f_conv_si( self, add=0.0 ):
        return conv_si_nu( self.f_order() + add )

    def f_order( self ):
        for_order   = sum( self.reactant_orders )
        for_order  += (1.0 if self.reaction_type == ReactionType.THREE_BODY else 0.0)
        return for_order

    def reaction_string_short( self ):
        delimiter   = '<=>' if self.is_reversible else '=>'
        return ''.join( [self._species_string( self.reactants ), delimiter, self._species_string( self.products )] )

    def reaction_string( self ):
        delimiter   = '<=>' if self.is_reversible else '=>'
        return ''.join( [self._stoich_species_string( self.reactants ), delimiter, self._stoich_species_string( self.products )] )

    def r_conv_si( self, add=0.0 ):
        return conv_si_nu( self.r_order() + add )

    def r_order( self ):
        rev_order   = sum( self.product_orders )
        rev_order  += (1.0 if self.reaction_type == ReactionType.THREE_BODY else 0.0)
        return rev_order

    def _add_species_string( self ):
        if self.reaction_type == ReactionType.DEFAULT:
            return ''
        elif self.reaction_type == ReactionType.THREE_BODY:
            return '+M'
        elif self.reaction_type == ReactionType.PRESSURE:
            return '(+{:})'.format( self.falloff_species )

    def _checklist( self, elements, specii, reaction_strings ):
        string          = self.reaction_string_short()
        n_duplicates    = reaction_strings.count( string )
        is_duplicate    = self.is_duplicate()
        last_troe       = (self.troe_coeff[3] != 0.0) if self.troe_coeff and len( self.troe_coeff ) == 4 else True

        return [
            (len( self.reactants.keys() ) > 0, 'No reactants are speciefied'),
            (len( self.products.keys() ) > 0, 'No products are speciefied'),
            (last_troe, 'Last TROE coefficient must not be zero (set to SMALL)'),
            (not (n_duplicates == 1 and is_duplicate), 'Reaction has duplicate flag, but no duplicate reaction is found'),
            (not (n_duplicates > 1 and not is_duplicate), 'Reaction is duplicate, but no duplicate flag is provided'),
            (check_element_balance( elements, specii, self ), 'Reaction does not balance in elements'),
            (check_mass_balance( specii, self ), 'Reaction does not balance in mass'),
        ]

    def _species_string( self, values ):
        tmp = '+'.join( [key for key in values] )
        return tmp + self._add_species_string()

    def _stoich_species_string( self, values ):
        tmp = '+'.join( ['{:}{:}'.format( as_short( values[key], round_digits=12 ), key ) for key in values] )
        return tmp + self._add_species_string()

    def __str__( self ):
        return self.reaction_string()

class ReactionContainer( BaseListContainer ):
    """Container (storage) object for Reaction class objects."""
    _type = Reaction

    def __init__( self, items=None, all_flag=False ):
        super().__init__( items=items, all_flag=all_flag )

        self.unit_k0    = def_unit_k0
        self.unit_Ea    = def_unit_Ea

    def chemkinify( self, unit_k0=def_unit_k0, unit_Ea=def_unit_Ea ):
        reaction_strings    = [x.reaction_string() for x in self]
        max_len             = max( len( x ) for x in reaction_strings ) + 1

        return [
            '! --- Reaction #{:} ---\n{:}\n'.format( idx+1,
                x.chemkinify( min_len=max_len, reaction_string=reaction_strings[idx], unit_k0=unit_k0, unit_Ea=unit_Ea )
            ) for idx, x in enumerate( self )
        ]

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
def cmsk_to_si( coeff, si_conv, unit_k0, unit_Ea ):
    """Convert the given coefficients to the specified units."""
    return [
        coeff[0] * si_conv * conv_k0[unit_k0],
        coeff[1],
        coeff[2] * conv_Ea[unit_Ea],
    ]

def si_to_cmsk( coeff, si_conv, unit_k0, unit_Ea ):
    """Convert the given coefficients to the specified units."""
    return [
        coeff[0] / si_conv / conv_k0[unit_k0],
        coeff[1],
        coeff[2] / conv_Ea[unit_Ea],
    ]

def conv_si_nu( sum_nu ):
    return CM_M ** (3 * (sum_nu - 1))

def check_element_balance( elements, specii, reaction, limit=1e-5 ):
    """Check the given reaction for their element balance."""
    elm_sum = {key: 0.0 for key in elements.keys()}

    for x, y in reaction.reactants.items():
        for key in specii[x].thermo.composition:
            elm_sum[key] += specii[x].thermo.composition[key] * y

    for x, y in reaction.products.items():
        for key in specii[x].thermo.composition:
            elm_sum[key] -= specii[x].thermo.composition[key] * y

    if any( abs( x ) > limit for x in elm_sum.values() ):
        logging.warning( 'Element conservation violation: {:}'.format( elm_sum ) )
    return all( abs( x ) <= limit for x in elm_sum.values() )

def check_mass_balance( specii, reaction, limit=1e-5 ):
    """Check the given reaction for their mass balance."""
    mass_sum = 0.0

    for x, y in reaction.reactants.items():
        mass_sum += specii[x].molar_mass * y

    for x, y in reaction.products.items():
        mass_sum -= specii[x].molar_mass * y

    if abs( mass_sum ) > limit:
        logging.warning( 'Mass conservation violation: {:}'.format( mass_sum ) )
    return abs( mass_sum ) <= limit