# Copyright 2025 NWChemEx Community
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
"""
@author: Felix Rojas
"""

import numpy as np
import pluginplay as pp
import tensorwrapper as tw
from simde import TotalEnergy


class LennardJonesPotential(pp.ModuleBase):
    # Module Construct --------------------------------------------------------
    def __init__(self):
        """
        This module Evaluates the Lennard-Jones 1D potential function (E)
        """
        pp.ModuleBase.__init__(self)
        self.description("Lennard-Jones 1D potential function")
        self.satisfies_property_type(TotalEnergy())

    # --------------------------------------------------------------------------

    # Module run_ member function ---------------------------------------------
    def run_(self, inputs, submods):
        """
        Parameters
        ----------
        inputs : Diatomic distance,
        TYPE ---> Float

        Returns
        -------
        E: Lennard-Jonnes 1D potential Energy,
        TYPE ---> Float
        """
        pt = TotalEnergy()
        (chem_sys,) = pt.unwrap_inputs(inputs)
        mol = chem_sys.molecule
        coor_0 = np.array([mol.at(0).x, mol.at(0).y, mol.at(0).z])
        coor_1 = np.array([mol.at(1).x, mol.at(1).y, mol.at(1).z])
        # ----------------------------------------------------------------------
        assert mol.size() == 2  # <--- To check molcule size contains 2-atoms
        # ----------------------------------------------------------------------
        r = np.linalg.norm(coor_0 - coor_1)
        # -------------- LENNARD-JONES FUNCTION --------------------------------
        E = 4 * ((1 / r**12) - (1 / r**6))
        # ------------- ANALYTIC FORCE -----------------------------------------
        # DE_x = -24 * ((2 / r**13) - (1 / r**7))
        # FC = -DE_x
        # ----------------------------------------------------------------------
        E = tw.Tensor(np.array(E))
        rv = self.results()
        return pt.wrap_results(rv, E)

    # --------------------------------------------------------------------------


def load_lennard_jones_potential(mm):
    mm.add_module("Lennard-Jones", LennardJonesPotential())
