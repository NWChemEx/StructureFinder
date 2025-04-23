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

import pluginplay as pp
from simde import (
    EnergyNuclearGradientStdVectorD,
    TotalEnergyNuclearOptimization,
    MoleculeFromString,
    TotalEnergy,
)
from berny import Berny, geomlib
import chemist
import numpy as np


class GeomoptViaPyberny(pp.ModuleBase):

    def __init__(self):
        pp.ModuleBase.__init__(self)
        self.satisfies_property_type(TotalEnergyNuclearOptimization())
        self.description("Performs PyBerny optimization")
        self.add_submodule(TotalEnergy(), "Energy")
        self.add_submodule(EnergyNuclearGradientStdVectorD(), "Gradient")
        self.add_submodule(MoleculeFromString(), "StringConv")

    def run_(self, inputs, submods):
        pt = TotalEnergyNuclearOptimization()
        sys, points = pt.unwrap_inputs(inputs)
        molecule = sys.molecule

        # Convert Chemist Chemical System to XYZ
        xyz = ""
        xyz += str(molecule.size()) + "\n\n"
        for i in range(molecule.size()):
            xyz += (molecule.at(i).name + " " + str(molecule.at(i).x) + " " +
                    str(molecule.at(i).y) + " " + str(molecule.at(i).z) + "\n")

        # Loads the geometry string into the Berny optimizer
        # object.
        optimizer = Berny(geomlib.loads(xyz, fmt="xyz"))

        for geom in optimizer:
            # Converts the "Berny" geometry object to Chemical System
            geom2xyz = geom.dumps("xyz")
            lines = geom2xyz.split("\n")
            mol_string = "\n".join(lines[2:])
            xyz2chem_mol = submods["StringConv"].run_as(
                MoleculeFromString(), mol_string)
            geom = chemist.ChemicalSystem(xyz2chem_mol)
            geom_nuclei = geom.molecule.nuclei.as_nuclei()
            geom_points = geom_nuclei.charges.point_set.as_point_set()

            # Main optimizer operation
            energy = submods["Energy"].run_as(TotalEnergy(), geom)
            gradients = submods["Gradient"].run_as(
                EnergyNuclearGradientStdVectorD(), geom, geom_points)
            optimizer.send((np.array(energy).item(), gradients))

        opt_geom_nuclei = geom.molecule.nuclei.as_nuclei()
        opt_geom_points = opt_geom_nuclei.charges.point_set.as_point_set()
        # Optimized energy is of type "float"
        e = energy
        print(e)
        rv = self.results()
        return pt.wrap_results(rv, e, opt_geom_points)


def load_pyberny_modules(mm):
    mm.add_module("PyBerny", GeomoptViaPyberny())
