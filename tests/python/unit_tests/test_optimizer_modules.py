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
from structurefinder.optimizer_modules import Scipy_optimizer_aoenergy, Scipy_optimizer_energy
import unittest
from nwchemex import load_modules
from pluginplay import ModuleManager
from simde import AOEnergy, Energy, MolecularBasisSet
from chemist import AOSpaceD, ChemicalSystem

class Test_scipy_optimizer_aoenergy(unittest.TestCase):
    def setUp(self):
        self.mol = ct.get_molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
        self.mm = ModuleManager()
        load_modules(self.mm)
        self.mod_key = 'scipy_optimizer_aoenergy'
        self.mm.add_module(self.mod_key, Scipy_optimizer_aoenergy())

    def test_scf(self):
        self.mm.change_submod(self.mod_key, 'AOEnergy', 'SCF Energy')
        bs = self.mm.at("sto-3g").run_as(MolecularBasisSet(), self.mol)
        aos = AOSpaceD(bs)
        chem_sys = ChemicalSystem(self.mol)
        energy = self.mm.at(self.mod_key).run_as(AOEnergy(), aos, chem_sys)
        self.assertAlmostEqual(energy, -1.0787343257869195, places=5)


class Test_scipy_optimizer_energy(unittest.TestCase):
    def setUp(self):
        self.mol = ct.get_molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
        self.mm = ModuleManager()
        load_modules(self.mm)
        self.mod_key = 'scipy_optimizer_energy'
        self.mm.add_module(self.mod_key, Scipy_optimizer_energy())

    def test_scf(self):
        self.mm.change_submod(self.mod_key, 'AOEnergy', 'SCF Energy')
        self.mm.change_submod(self.mod_key, 'MolecularBasisSet', 'sto-3g')
        chem_sys = ChemicalSystem(self.mol)
        energy = self.mm.at(self.mod_key).run_as(Energy(), chem_sys)
        self.assertAlmostEqual(energy, -1.0787343257869195, places=5)
