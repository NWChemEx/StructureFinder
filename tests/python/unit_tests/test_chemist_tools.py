# Copyright 2023 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import structurefinder.chemist_tools as ct
import unittest

class Test_chemist_tools(unittest.TestCase):    
    def setUp(self):
        self.mol = ct.get_molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])

    def test_get_atomic_mass(self):
        self.assertEqual(ct.get_atomic_mass(1), 1.0079)
        self.assertEqual(ct.get_atomic_mass('H'), 1.0079)
        self.assertEqual(ct.get_atomic_mass(6), 12.011)
        self.assertEqual(ct.get_atomic_mass('C'), 12.011)

    def test_get_atomic_number(self):
        self.assertEqual(ct.get_atomic_number('H'), 1)
        self.assertEqual(ct.get_atomic_number('C'), 6)

    def test_get_molecule(self):
        self.assertEqual(self.mol.size(), 2)
        self.assertEqual(self.mol.at(0).name, 'H')
        self.assertEqual(self.mol.at(1).name, 'H')
        self.assertEqual(self.mol.at(0).x, 0.0)
        self.assertEqual(self.mol.at(0).y, 0.0)
        self.assertEqual(self.mol.at(0).z, 0.0)
        self.assertEqual(self.mol.at(1).x, 0.0)
        self.assertEqual(self.mol.at(1).y, 0.0)
        self.assertEqual(self.mol.at(1).z, 0.74) 
 
    def test_get_molecule_coordinates(self):
        coords = ct.get_molecule_coordinates(self.mol)
        self.assertEqual(coords, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])

    def test_get_molecule_symbols(self):
        symbols = ct.get_molecule_symbols(self.mol)
        self.assertEqual(symbols, ['H', 'H'])

    def test_get_periodic_table(self):
        pt = ct.get_periodic_table()
        self.assertEqual(len(pt), 54)
        self.assertEqual(pt[1], 'H')
        self.assertEqual(pt[6], 'C')
    
    def test_get_symbol(self):
        self.assertEqual(ct.get_symbol(1), 'H')
        self.assertEqual(ct.get_symbol(6), 'C')