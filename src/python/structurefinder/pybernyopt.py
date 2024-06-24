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

from berny import Berny, geomlib


def optimize_pyberny(molecule):
    """
    This function takes a chemical system and uses PyBerny
    to optimize the provided geometry (molecule) and
    returns the energy of the optimized system.
    """
    # Convert a Chemical System to an XYZ coordinate string
    xyz = ""
    xyz += (str(molecule.size()) + "\n\n")
    for i in range(molecule.size()):
        xyz += (molecule.at(i).name + " " + str(molecule.at(i).x) + " " +
                str(molecule.at(i).y) + " " + str(molecule.at(i).z) + "\n")

    # Loads the geometry string into the Berny optimizer
    # object.
    optimizer = Berny(geomlib.loads(xyz, fmt='xyz'))

    for geom in optimizer:
        energy = calculate_energy(geom)
        gradients = calculate_gradient(geom)
        optimizer.send((energy, gradients))

    relaxed = geom
    xyz_opt = relaxed.dumps(fmt='xyz')

    # Optimized energy is of type "float"
    return energy, xyz_opt
