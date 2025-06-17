# # Copyright 2024 NWChemEx Community
# #
# # Licensed under the Apache License, Version 2.0 (the "License");
# # you may not use this file except in compliance with the License.
# # You may obtain a copy of the License at
# #
# # http://www.apache.org/licenses/LICENSE-2.0
# #
# # Unless required by applicable law or agreed to in writing, software
# # distributed under the License is distributed on an "AS IS" BASIS,
# # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# # See the License for the specific language governing permissions and
# # limitations under the License.

import structurefinder
import nwchemex
import pluginplay as pp
import chemist
import unittest
from simde import TotalEnergyNuclearOptimization
import numpy as np
import tensorwrapper as tw

def print_pointset(pointset):
    printout = ' '
    for i in range(pointset.size()):
        printout+= '['
        for j in range(3):
            printout+= str(pointset.at(i).coord(j)) + ' '
        printout+= ']'
    print(printout)




class Test_TotalEnergyNuclearOptimization(unittest.TestCase):

    def test_optimize_BEfire(self):
        mm = pp.ModuleManager()
        nwchemex.load_modules(mm)
        structurefinder.load_modules(mm)
        mm.change_input("NWChem : SCF", "basis set", "sto-3g")
        mm.change_input("NWChem : SCF Gradient", "basis set", "sto-3g")
        mm.change_submod("BackwardEulerFire", "Gradient", "NWChem : SCF Gradient")
        mm.change_submod("BackwardEulerFire", "Energy", "NWChem : SCF")
        mm.change_submod("BackwardEulerFire", "StringConv",
                         "ChemicalSystem via QCElemental")
        egy, pts = mm.run_as(TotalEnergyNuclearOptimization(), "BackwardEulerFire",
                       self.sys, self.pointset)
        print("Energy = " + str(egy))
        print_pointset(pts)
        print(self.sys.molecule)
        # print(pts) <-- chemist types missing python string representation
        #self.assertAlmostEqual(np.array(egy), -1.117505879316, 10)

    def setUp(self):
        self.mol = chemist.Molecule()
        self.mol.push_back(chemist.Atom("H", 1, 1.0079, 0.0, 0.0, 0.0))
        self.mol.push_back(chemist.Atom("H", 1, 1.0079, 0.0, 0.0, 1.0))
        self.sys = chemist.ChemicalSystem(self.mol)
        self.nuclei = self.mol.nuclei.as_nuclei()
        self.pointset = self.nuclei.charges.point_set.as_point_set()
