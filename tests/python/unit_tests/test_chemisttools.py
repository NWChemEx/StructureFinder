import structurefinder.chemisttools as ct
import unittest

class Testpt2driver(unittest.TestCase):

    def test_get_coordinates(self):
        mol = ct.get_molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
        coords = ct.get_coordinates(mol)
        self.assertEqual(coords, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])

    def test_get_symbols(self):
        mol = ct.get_molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
        symbols = ct.get_symbols(mol)
        self.assertEqual(symbols, ['H', 'H'])

    def test_get_molecule(self):
        mol = ct.get_molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
        self.assertEqual(mol.size(), 2)
        self.assertEqual(mol.at(0).name, 'H')
        self.assertEqual(mol.at(1).name, 'H')
        self.assertEqual(mol.at(0).x, 0.0)
        self.assertEqual(mol.at(0).y, 0.0)
        self.assertEqual(mol.at(0).z, 0.0)
        self.assertEqual(mol.at(1).x, 0.0)
        self.assertEqual(mol.at(1).y, 0.0)
        self.assertEqual(mol.at(1).z, 0.74)        

    def test_get_symbol(self):
        self.assertEqual(ct.get_symbol(1), 'H')

    def test_get_atomic_number(self):
        self.assertEqual(ct.get_atomic_number('H'), 1)

