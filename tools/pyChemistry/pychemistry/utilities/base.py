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
    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __str__( self ):
        raise( NotImplementedError )

    def _checklist( self, **kwargs ):
        return [(False, 'No checks defined')]

    def is_valid( self, **kwargs ):
        warnings = [warning for check, warning in self._checklist( **kwargs ) if not check]

        # provide logging warnings for each failed check
        for warning in warnings:
            logging.warning( '{:}:{:}'.format( self, warning ) )

        return (not warnings)

class BaseOrderedDictContainer( object ):
    """Container (storage) object for Base class objects."""
    _type = Base

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
            raise(TypeError(BaseOrderedDictContainer, type(other)))

    def __getitem__(self, key):
        return self._items.__getitem__(key)

    def __init__( self, items=None, all_flag=False ):
        if items is not None:
            for item in items.values():
                self._validate_type(item)

        self._items = OrderedDict() if items is None else items
        self.all_flag = all_flag

    def __iter__(self):
        return self._items.__iter__()

    def __len__(self):
        return self._items.__len__()

    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __setitem__(self, key, value):
        self._validate_type(value)
        return self._items.__setitem__(key, value)

    def __str__( self ):
        return str( self.__dict__ )

    def _combine( self, other ):
        for key, value in other.items():
            self._validate_type( value )

            if key not in self:
                self[key] = value
            else:
                logging.debug( 'Ignore duplicate data for "{:}"!'.format( key ) )

    def _validate_type( self, item ):
        if self._type != type( item ):
            raise(TypeError(self._type, type(item)))

    def chemkinify( self, keys=None, **kwargs ):
        tmp = self.keys() if keys is None else keys
        return ['{:}\n'.format( self[key].chemkinify( **kwargs ) ) for key in tmp]

    def is_valid( self, keys=None, **kwargs ):
        tmp = self.keys() if keys is None else keys
        return all( [self[key].is_valid( **kwargs ) for key in tmp] )

    def items(self):
        return self._items.items()

    def keys(self):
        return self._items.keys()

    def values(self):
        return self._items.values()

class BaseListContainer( object ):
    """Container (storage) object for Base class objects."""
    _type = Base

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
            raise(TypeError(BaseListContainer, type(other)))

    def __init__( self, items=None, all_flag=False ):
        if items is not None:
            for item in items:
                self._validate_type(item)

        self._items = [] if items is None else items
        self.all_flag = all_flag

    def __iter__(self):
        return self._items.__iter__()

    def __len__(self):
        return self._items.__len__()

    def __repr__( self ):
        return '<{:}: {:}>'.format( self.__class__.__name__, self )

    def __str__( self ):
        return str( self.__dict__ )

    def _combine( self, other ):
        for value in other:
            self._validate_type( value )
            self.append( other )

    def _validate_type( self, item ):
        if self._type != type( item ):
            raise(TypeError(self._type, type(item)))

    def append(self, item):
        self._validate_type(item)
        self._items.append(item)

    def chemkinify( self, **kwargs ):
        return ['{:}\n'.format( x.chemkinify( **kwargs ) ) for x in self]

    def is_valid( self, **kwargs ):
        return all( [x.is_valid( **kwargs ) for x in self] )

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------