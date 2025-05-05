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

import structurefinder
import nwchemex
import numpy as np
import pluginplay as pp
import chemist
import unittest
from simde import TotalEnergyNuclearOptimization


def diatomic_bond_distance(coords):
    val = 0
    for i in range(int(len(coords) / 2)):
        val += (coords[i] - coords[i + 3])**2
    distance = np.sqrt(val)
    return distance


class Test_optimize_pyberny(unittest.TestCase):

    def test_optimize_pyberny(self):
        mm = pp.ModuleManager()
        nwchemex.load_modules(mm)
        structurefinder.load_modules(mm)

        pyberny_mod = mm.at("PyBerny")
        nwchem_scf_mod = mm.at("NWChem : SCF")
        nwchem_grad_mod = mm.at("NWChem : SCF Gradient")
        string_conv_mod = mm.at("ChemicalSystem via QCElemental")

        nwchem_scf_mod.change_input("basis set", "sto-3g")
        nwchem_grad_mod.change_input("basis set", "sto-3g")
        pyberny_mod.change_submod("Energy", nwchem_scf_mod)
        pyberny_mod.change_submod("Gradient", nwchem_grad_mod)
        pyberny_mod.change_submod("StringConv", string_conv_mod)

        energy, points = pyberny_mod.run_as(TotalEnergyNuclearOptimization(),
                                            self.sys, self.point_set_i)

        coords = []
        for atom in range(points.size()):
            for coord in range(3):
                coords.append(points.at(atom).coord(coord))

        distance = diatomic_bond_distance(coords)

        energy = np.array(energy).item()
        self.assertAlmostEqual(energy, self.energy, 10)
        self.assertAlmostEqual(distance, self.distance, 8)

    def setUp(self):
        self.mol = chemist.Molecule()
        self.mol.push_back(chemist.Atom("H", 1, 1.0079, 0.0, 0.0, 0.0))
        self.mol.push_back(chemist.Atom("H", 1, 1.0079, 0.0, 0.0, 1.0))
        self.sys = chemist.ChemicalSystem(self.mol)
        self.nuclei = self.mol.nuclei.as_nuclei()
        self.point_set_i = self.nuclei.charges.point_set.as_point_set()
        self.distance = 1.34606231
        self.energy = -1.117505879316
