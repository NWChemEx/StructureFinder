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
from simde import EnergyNuclearGradientStdVectorD, TotalEnergy, MoleculeFromString
from berny import Berny, geomlib, optimize
import chemist
import qcelemental as qcel


class GeomoptViaPyberny(pp.ModuleBase):

    def __init__(self):
        pp.ModuleBase.__init__(self)
        self.satisfies_property_type(TotalEnergy())
        self.description("Performs PyBerny optimization")
        self.add_submodule(TotalEnergy(), "Energy")
        self.add_submodule(EnergyNuclearGradientStdVectorD(), "Gradient")
        self.add_submodule(MoleculeFromString(), "StringConv")
        
    def run_(self, inputs, submods):
        pt = TotalEnergy()
        mol, = pt.unwrap_inputs(inputs)
        molecule = mol.molecule

        # Convert Chemist Chemical System to XYZ
        xyz = ""
        xyz += (str(molecule.size()) + "\n\n")
        for i in range(molecule.size()):
            xyz += (molecule.at(i).name + " " + str(molecule.at(i).x) + " " +
                    str(molecule.at(i).y) + " " + str(molecule.at(i).z) + "\n")

        # Loads the geometry string into the Berny optimizer
        # object.
        optimizer = Berny(geomlib.loads(xyz, fmt='xyz'))

        for geom in optimizer:

            # Converts the "Berny" geometry object to Chemical System
            geom2xyz = geom.dumps('xyz')
            print('Berny Geom to XYZ value: \n' + geom2xyz + '\n')
            lines = geom2xyz.split('\n')
            print('Lines of geom2xyz: \n' + str(lines) + '\n')
            mol_string = '\n'.join(lines[2:])
            print('Lines to string: \n' + mol_string + '\n')
            xyz2chem_mol = submods["StringConv"].run_as(MoleculeFromString(), mol_string)
            print('String conversion from xyz to chem sys: \n' + str(xyz2chem_mol.nuclei) + '\n')
            geom = chemist.ChemicalSystem(xyz2chem_mol)
            print('Chemical system of xyz2chem_mol: \n' + str(geom.molecule.nuclei) + '\n')

            # Main optimizer operation
            energy = submods["Energy"].run_as(TotalEnergy(), geom)
            print('Interim energy: \n' + str(energy) + '\n')
            gradients = submods["Gradient"].run_as(
                EnergyNuclearGradientStdVectorD(), geom, chemist.PointSetD())
            print('Interim gradient: \n' + str(gradients) + '\n')
            optimizer.send((energy, gradients))

        opt_geom = geom.molecule.nuclei
        print('Resulting relaxed geometry (assigned to variable opt_geom): \n' + str(opt_geom))
        # Optimized energy is of type "float"
        e = energy
        print(e)
        rv = self.results()
        return pt.wrap_results(rv, e)


def load_pyberny_modules(mm):
    mm.add_module("PyBerny", GeomoptViaPyberny())
