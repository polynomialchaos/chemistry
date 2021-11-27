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


class Base():
    """Base object (storing data)."""

    def __repr__(self):
        return '<{:}: {:}>'.format(self.__class__.__name__, self)

    def __str__(self):
        return str(self.__dict__)

    def _check_list(self, **_):
        return [(False, 'Missing checks for "{:}" class!'.format(
            self.__class__.__name__
        ))]

    def chemkinify(self):
        """Return a CHEMKIN formatted string."""
        raise NotImplementedError('chemkinify')

    def is_valid(self, **kwargs):
        """Check if the object is valid (runs check list)."""
        errors = [w for c, w in self._check_list(**kwargs) if not c]

        # provide logging warnings for each failed check
        for err in errors:
            logging.warning('%s:%s', self, err)

        return not errors


class _BaseContainer():
    """Base container object (storing datas)."""
    _store_type = Base

    def __add__(self, other):
        if isinstance(other, type(self)):
            if other is None:
                return deepcopy(self)

            if self.all_flag or other.all_flag:
                raise Exception(
                    'Addition of data with ALL flag set is not supported!')

            new = deepcopy(self)
            new._combine(other)
            return new

        raise TypeError('Wrong type "{:}" ({:})!'.format(
            type(other), type(self)))

    def __init__(self, all_flag=False):
        self.all_flag = all_flag

    def __repr__(self):
        return '<{:}: {:}>'.format(self.__class__.__name__, self)

    def __str__(self):
        return str(self.__dict__)

    def _combine(self, other):
        raise NotImplementedError('_combine')

    def _validate_item(self, item):
        if self._store_type != type(item):
            raise TypeError('Wrong type "{:}" ({:})!'.format(
                type(item), type(self._store_type)))


class BaseDictContainer(_BaseContainer):
    """Base dict container object (storing datas)."""

    def __getitem__(self, key):
        return self._items.__getitem__(key)

    def __init__(self, items=None, all_flag=False):
        super().__init__(all_flag=all_flag)
        self._items = {}

        if items is not None:
            for key, value in items.items():
                self[key] = value

    def __iter__(self):
        return self._items.__iter__()

    def __len__(self):
        return self._items.__len__()

    def __setitem__(self, key, value):
        self._validate_item(value)
        return self._items.__setitem__(key, value)

    def _combine(self, other):
        for key, value in other.items():
            if key in self:
                logging.debug('Ignore duplicate data for "%s"!', key)
                continue

            self[key] = value

    def chemkinify(self, keys=None, **kwargs):
        """Return a list of CHEMKIN formatted strings."""
        tmp = self.keys() if keys is None else keys
        return ['{:}\n'.format(self[key].chemkinify(**kwargs)) for key in tmp]

    def is_valid(self, keys=None, **kwargs):
        """Check if the objects in container ar valid."""
        tmp = self.keys() if keys is None else keys
        return all([self[key].is_valid(**kwargs) for key in tmp])

    def items(self):
        """Return dictionary items."""
        return self._items.items()

    def keys(self):
        """Return dictionary keys."""
        return self._items.keys()

    def values(self):
        """Return dictionary values."""
        return self._items.values()


class BaseListContainer(_BaseContainer):
    """Base list container object (storing datas)."""

    def __init__(self, items=None, all_flag=False):
        super().__init__(all_flag=all_flag)
        self._items = []

        if items is not None:
            for item in items:
                self.append(item)

    def __iter__(self):
        return self._items.__iter__()

    def __len__(self):
        return self._items.__len__()

    def _combine(self, other):
        for value in other:
            self.append(value)

    def append(self, item):
        """Append a item to the list."""
        self._validate_item(item)
        self._items.append(item)

    def chemkinify(self, **kwargs):
        """Return a list of CHEMKIN formatted strings."""
        return ['{:}\n'.format(x.chemkinify(**kwargs)) for x in self]

    def is_valid(self, **kwargs):
        """Check if the objects in container ar valid."""
        return all([x.is_valid(**kwargs) for x in self])
