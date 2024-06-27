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
from simde import TotalEnergy, EnergyNuclearGradientStdVectorD
from berny import Berny, geomlib, optimize


class GeomoptViaPyberny(pp.ModuleBase):

    def __init__(self):
        pp.ModuleBase.__init__(self)
        self.satisfies_property_type(EnergyNuclearGradientStdVectorD())
        self.description("Performs PyBerny optimization")
        self.add_submodule(EnergyNuclearGradientStdVectorD(), "Energy and Gradient")

    def run_(self, inputs, submods):
        pt = EnergyNuclearGradientsStdVectorD()
        mol, pointset1 = pt.unwrap_inputs(inputs)
        molecule = mol.molecule

        xyz = ""
        xyz += (str(molecule.size()) + "\n\n")
        for i in range(molecule.size()):
            xyz += (molecule.at(i).name + " " + str(molecule.at(i).x) + " " +
                    str(molecule.at(i).y) + " " + str(molecule.at(i).z) + "\n")

        # Loads the geometry string into the Berny optimizer
        # object.
        optimizer = Berny(geomlib.loads(xyz, fmt='xyz'))

        for geom in optimizer:
            xyz2qc_mol = qcel.models.Molecule.from_data(geom.dumps('xyz'))
            qc_mol2chemicalsystem = chemical_system_conversions.qc_mol2molecule(
                xyz2qc_mol)
            geom = chemist.ChemicalSystem(qc_mol2chemicalsystem)
            energy, gradients = submods["Energy and Gradient"].run_as(
                EnergyNuclearGradientStdVectorD(), geom)
            optimizer.send((energy, gradients))

        relaxed = geom
        xyz_opt = relaxed.dumps(fmt='xyz')
        print(xyz_opt)
        # Optimized energy is of type "float"
        e = energy
        rv = self.results()
        return pr.wrap_results(rv, e)


def load_pyberny_modules(mm):
    mm.add_module("PyBerny", GeomoptViaPyberny())
