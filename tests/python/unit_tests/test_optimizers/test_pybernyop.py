# Copyright 2024 NWChemEx Community
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

from structurefinder import load_modules
import pluginplay as pp
import chemist
import unittest


class Test_optimize_pyberny(unittest.TestCase):

    def test_optimize_pyberny(self):
        mm = pp.ModuleManager()
        load_modules(mm)
        egy = mm.run_as("PyBerny", chemist.ChemicalSystem(self.mol))
        print("Energy = " + egy)

    def setUp(self):
        self.mol = chemist.Molecule()
        self.mol.push_back(chemist.Atom("H", 1, 1.0079, 0.0, 0.0, 0.0))
        self.mol.push_back(chemist.Atom("H", 1, 1.0079, 0.0, 0.0, 1.0))
