####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) Florian Eigentler | 2020
####################################################################################################################################
import logging
from collections import OrderedDict
from copy import deepcopy

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
class Base( object ):
    """Object for storing data."""
    def is_valid( self, **kwargs ):
        warnings = [warning for check, warning in self._checklist( **kwargs ) if not check]

        # provide logging warnings for each failed check
        for warning in warnings:
            logging.warning( '{:}:{:}'.format( self, warning ) )

        return (not warnings)

    def _checklist( self, **kwargs ):
        return [(False, 'No checks defined')]

    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __str__( self ):
        raise( NotImplementedError )

class BaseOrderedDictContainer( OrderedDict ):
    """Container (storage) object for Base class objects."""
    _type = Base

    def __init__( self, items=None, all_flag=False ):
        super().__init__( [] if items is None else items )

        self.all_flag = all_flag

    def chemkinify( self, keys=None, **kwargs ):
        tmp = self.keys() if keys is None else keys
        return ['{:}\n'.format( self.__getitem__( key ).chemkinify( **kwargs ) ) for key in tmp]

    def is_valid( self, keys=None, **kwargs ):
        tmp = self.keys() if keys is None else keys
        return all( [self.__getitem__( key ).is_valid( **kwargs ) for key in tmp] )

    def _combine( self, other ):
        for key, value in other.items():
            self._validate_type( value )

            if key not in self.keys():
                self.__setitem__( key, value )
            else:
                logging.debug( 'Ignore duplicate data for "{:}"!'.format( key ) )

    def _validate_type( self, item ):
        if self._type != type( item ):
            raise( TypeError( 'Provided item has wrong type "{:}" ("{:}")!'.format( type( item ), self._type ) ) )

    def __add__( self, other ):
        if isinstance( self, other.__class__ ):
            if other is None: return deepcopy( self )
            elif not self.keys(): return deepcopy( other )

            if self.all_flag or other.all_flag:
                raise( Exception( 'Addition of data with ALL flag set is not supported!' ) )
            else:
                new = deepcopy( self )
                new._combine( other )
                return new
        else:
            raise( TypeError( 'Provided container has wrong type "{:}" ("{:}")!'.format( self.__class__, other.__class__ ) ) )

    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __str__( self ):
        return str( self.__dict__ )

class BaseListContainer( list ):
    """Container (storage) object for Base class objects."""
    _type = Base

    def __init__( self, items=None, all_flag=False ):
        super().__init__( [] if items is None else items )

        self.all_flag = all_flag

    def chemkinify( self, **kwargs ):
        return ['{:}\n'.format( x.chemkinify( **kwargs ) ) for x in self]

    def is_valid( self, **kwargs ):
        return all( [x.is_valid( **kwargs ) for x in self] )

    def _combine( self, other ):
        for value in other:
            self._validate_type( value )
            self.append( other )

    def _validate_type( self, item ):
        if self._type != type( item ):
            raise( TypeError( 'Provided item has wrong type "{:}" ("{:}")!'.format( type( item ), self._type ) ) )

    def __add__( self, other ):
        if isinstance( self, other.__class__ ):
            if other is None: return deepcopy( self )
            elif not len( self ): return deepcopy( other )

            if self.all_flag or other.all_flag:
                raise( Exception( 'Addition of data with ALL flag set is not supported!' ) )
            else:
                new = deepcopy( self )
                new._combine( other )
                return new
        else:
            raise( TypeError( 'Provided container has wrong type "{:}" ("{:}")!'.format( self.__class__, other.__class__ ) ) )

    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __str__( self ):
        return str( self.__dict__ )

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------