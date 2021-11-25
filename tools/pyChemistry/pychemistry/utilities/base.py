################################################################################
# @file base.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import logging
from copy import deepcopy


class Base(object):
    """Base object (storing data)."""

    def __repr__(self):
        return '<{:}: {:}>'.format(self.__class__.__name__, self)

    def __str__(self):
        return str(self.__dict__)

    def _data_check(self, **kwargs):
        return [(False, 'No data check defined')]

    def is_valid(self, **kwargs):
        errors = [w for c, w in self._data_check(**kwargs) if not c]

        # provide logging warnings for each failed check
        for err in errors:
            logging.warning('{:}:{:}'.format(self, err))

        return (not errors)


class _BaseContainer(object):
    """Base container object (storing datas)."""
    _store_type = Base

    def __add__(self, other):
        if type(self) == type(other):
            if other is None:
                return deepcopy(self)

            if self.all_flag or other.all_flag:
                raise(Exception(
                    'Addition of data with ALL flag set is not supported!'))
            else:
                new = deepcopy(self)
                new._combine(other)
                return new
        else:
            raise(TypeError(type(self), type(other)))

    def __iter__(self):
        return self._items.__iter__()

    def __len__(self):
        return self._items.__len__()

    def __repr__(self):
        return '<{:}: {:}>'.format(self.__class__.__name__, self)

    def __str__(self):
        return str(self.__dict__)

    def _combine(self, other):
        raise(NotImplementedError('_combine'))

    def _validate_item(self, item):
        if self._store_type != type(item):
            raise(TypeError(self._store_type, type(item)))

    def chemkinify(self, **kwargs):
        raise(NotImplementedError('chemkinify'))

    def is_valid(self, **kwargs):
        raise(NotImplementedError('is_valid'))


class BaseDictContainer(_BaseContainer):
    """Base dict container object (storing datas)."""

    def __getitem__(self, key):
        return self._items.__getitem__(key)

    def __init__(self, items=None, all_flag=False):
        self._items = {}
        self.all_flag = all_flag

        if items is not None:
            for key, value in items.items():
                self[key] = value

    def __setitem__(self, key, value):
        self._validate_item(value)
        return self._items.__setitem__(key, value)

    def _combine(self, other):
        for key, value in other.items():
            if key in self:
                logging.debug('Ignore duplicate data for "{:}"!'.format(key))
                continue

            self[key] = value

    def chemkinify(self, keys=None, **kwargs):
        tmp = self.keys() if keys is None else keys
        return ['{:}\n'.format(self[key].chemkinify(**kwargs)) for key in tmp]

    def is_valid(self, keys=None, **kwargs):
        tmp = self.keys() if keys is None else keys
        return all([self[key].is_valid(**kwargs) for key in tmp])

    def items(self):
        return self._items.items()

    def keys(self):
        return self._items.keys()

    def values(self):
        return self._items.values()


class BaseListContainer(_BaseContainer):
    """Base list container object (storing datas)."""

    def __init__(self, items=None, all_flag=False):
        self._items = []
        self.all_flag = all_flag

        if items is not None:
            for item in items:
                self.append(item)

    def _combine(self, other):
        for value in other:
            self.append(value)

    def append(self, item):
        self._validate_item(item)
        self._items.append(item)

    def chemkinify(self, **kwargs):
        return ['{:}\n'.format(x.chemkinify(**kwargs)) for x in self]

    def is_valid(self, **kwargs):
        return all([x.is_valid(**kwargs) for x in self])
