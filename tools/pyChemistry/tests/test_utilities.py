################################################################################
# @file test_utilities.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
import unittest
from pychemistry.utilities import Base, BaseDictContainer, BaseListContainer, transport
from pychemistry.utilities import Element, ElementContainer, REF_ELEMENTS
from pychemistry.utilities import Species, SpeciesContainer, NA
from pychemistry.utilities import Thermo, ThermoContainer
from pychemistry.utilities import Transport, TransportContainer
from pychemistry.utilities import Reaction, ReactionContainer
from pychemistry.utilities import as_int, as_short, chunk_list, is_number
from pychemistry.utilities import species
from pychemistry.utilities.reaction import ReactionType
from pychemistry.utilities.species import SpeciesContainer


class TestBase(unittest.TestCase):

    def test_base(self):
        base = Base()
        self.assertEqual(base.is_valid(), False)

    def test_base_dict_container(self):
        elements = [('E1', Base()), ('E2', Base())]
        container = BaseDictContainer()
        self.assertRaises(TypeError, BaseDictContainer, {'key': 'value'})
        self.assertRaises(TypeError, container.__add__, 12)
        for k, v in elements:
            container[k] = v
        self.assertEqual(len(container), 2)
        self.assertListEqual(list(container.keys()), [x[0] for x in elements])
        self.assertListEqual(list(container.values()),
                             [x[1] for x in elements])
        self.assertListEqual(list(container.items()), [
                             (x[0], x[1]) for x in elements])

    def test_base_list_container(self):
        elements = [Base(), Base()]
        container = BaseListContainer()
        self.assertRaises(TypeError, BaseListContainer, [12])
        self.assertRaises(TypeError, container.__add__, 12)
        container.append(elements[0])
        container.append(elements[1])
        self.assertEqual(len(container), 2)


class TestElement(unittest.TestCase):
    def test_element(self):
        element = Element('H')
        self.assertEqual(element.mass, REF_ELEMENTS['H'])
        element.mass = 13.0
        self.assertEqual(element.mass, 13.0)
        element.chemkinify()

    def test_element_container(self):
        elements = [Element('H'), Element('H2')]
        container = ElementContainer({elements[0].symbol: elements[0]})
        for element in elements:
            container[element.symbol] = element
        self.assertEqual(len(container), 2)
        container.chemkinify()


class TestReaction(unittest.TestCase):
    def test_reaction(self):
        reaction = Reaction(ReactionType.DEFAULT, reactants={'H2O': 1},
                            products={'H2O': 1}, is_reversible=True, arr_coeff=(1,2,3))
        reaction.chemkinify()

    def test_reaction_container(self):
        elements = [Reaction(ReactionType.DEFAULT, reactants={'H2O': 1},
                             products={'H2O': 1}, is_reversible=True, arr_coeff=(1,2,3)),
                    Reaction(ReactionType.DEFAULT, reactants={'H2O': 1},
                             products={'H2O': 1}, is_reversible=True, arr_coeff=(1,2,3)),
                    ]
        container = ReactionContainer(elements)
        self.assertEqual(len(container), 2)


class TestSpecies(unittest.TestCase):
    def test_species(self):
        species = Species('H')
        self.assertRaises(ValueError, getattr, species, 'molecule_weight')
        self.assertRaises(ValueError, getattr, species, 'molar_mass')
        self.assertRaises(AttributeError, setattr,
                          species, 'molecule_weight', 12.0)
        species.molar_mass = 12.0
        self.assertEqual(species.molar_mass, 12.0)
        self.assertEqual(species.molecule_weight, 12.0 / NA)
        species.chemkinify()

    def test_species_container(self):
        elements = [Species('H'), Species('H2')]
        container = SpeciesContainer({elements[0].symbol: elements[0]})
        for element in elements:
            container[element.symbol] = element
        self.assertEqual(len(container), 2)
        container.chemkinify()


class TestThermo(unittest.TestCase):
    def test_thermo(self):
        thermo = Thermo('H', info='test', composition={'H': 1},
                        phase='G', bounds=(1, 2, 3), coeff_low=(1, 2, 3, 4, 5, 6, 7),
                        coeff_high=(8, 9, 10, 11, 12, 13, 14))
        thermo.chemkinify()

    def test_thermo_container(self):
        elements = [Thermo('H'), Thermo('H2')]
        container = ThermoContainer({elements[0].symbol: elements[0]})
        for element in elements:
            container[element.symbol] = element
        self.assertEqual(len(container), 2)


class TestTransport(unittest.TestCase):
    def test_transport(self):
        transport = Transport('H', geom=0, pot_lj=1,
                              col_lj=1, dip_mo=1, pol=1, rot_rel=1)
        transport.chemkinify()

    def test_transport_container(self):
        elements = [Transport('H'), Transport('H2')]
        container = TransportContainer({elements[0].symbol: elements[0]})
        for element in elements:
            container[element.symbol] = element
        self.assertEqual(len(container), 2)


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
