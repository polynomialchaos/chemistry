################################################################################
# @file test_utilities.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import unittest
from pychemistry.utilities import Base, BaseDictContainer, BaseListContainer
from pychemistry.utilities import as_int, as_short, chunk_list, is_number


class TestBase(unittest.TestCase):

    def test_base_object(self):
        base = Base()
        self.assertEqual(base.is_valid(), False)

    def test_base_dict_container(self):
        base_container = BaseDictContainer()
        self.assertRaises(TypeError, BaseDictContainer, {'key': 'value'})
        self.assertRaises(TypeError, base_container.__add__, 12)
        base_container['test'] = Base()

    def test_base_list_container(self):
        base_container = BaseListContainer()
        self.assertRaises(TypeError, BaseListContainer, [12])
        self.assertRaises(TypeError, base_container.__add__, 12)
        base_container.append(Base())


class TestUtilities(unittest.TestCase):

    def test_as_int(self):
        self.assertEqual(as_int(0.0), 0)
        self.assertEqual(as_int(1.0), 1)
        self.assertNotEqual(as_int(0.1), 0)
        self.assertNotEqual(as_int(1.1), 1)
        self.assertEqual(as_int(0.1, round_digits=0), 0)
        self.assertEqual(as_int(1.1, round_digits=0), 1)

    def test_as_short(self):
        self.assertEqual(as_short(1), '')
        self.assertEqual(as_short(0.0), '0')
        self.assertEqual(as_short(1.1), '1.1')
        self.assertEqual(as_short(1.1, round_digits=0), '')

    def test_chunk_list(self):
        elements = [1, 2, 3, 4, 5]
        self.assertListEqual([x for x in chunk_list(elements, 1)], [
                             [1], [2], [3], [4], [5]])
        self.assertListEqual([x for x in chunk_list(elements, 2)], [
                             [1, 2], [3, 4], [5]])
        self.assertListEqual([x for x in chunk_list(
            elements, len(elements))], [elements])
        self.assertListEqual([x for x in chunk_list(
            elements, len(elements) + 1)], [elements])

    def test_is_number(self):
        self.assertEqual(is_number(1), True)
        self.assertEqual(is_number(1.0), True)
        self.assertEqual(is_number('1.0e16'), True)

        self.assertEqual(is_number(True), False)
        self.assertEqual(is_number('Hello World!'), False)


################################################################################
# CALL BY SCRIPT
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    unittest.main()
