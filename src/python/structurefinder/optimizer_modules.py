# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

from pluginplay import ModuleBase
from simde import AOEnergy, Energy, MolecularBasisSet
from chemist import AOSpaceD, ChemicalSystem
from . import chemist_tools as ct
from numpy import ravel
from scipy.optimize import minimize

class Scipy_optimizer_aoenergy(ModuleBase):
    """Optimization module based on scipy.optimize.minimize. 
    Satisfies simde.AOEnergy property type."""
    ptype = AOEnergy()

    def __init__(self):
        ModuleBase.__init__(self)
        self.description(self.__doc__)
        self.satisfies_property_type(self.ptype)
        self.add_submodule(AOEnergy(), "AOEnergy")

    def run_(self, inputs, submods):
        [aos, system_in] = self.ptype.unwrap_inputs(inputs)
        mol = system_in.molecule
        coords = ct.get_molecule_coordinates(mol)
        symbols = ct.get_molecule_symbols(mol)
        def _get_energy(coords):
            m = ct.get_molecule(symbols, coords)
            chem_sys = ChemicalSystem(m)
            E = submods["AOEnergy"].run_as(self.ptype, aos, chem_sys)
            return E
        result = minimize(_get_energy, ravel(coords))
        E = result.fun
        r = self.results()
        return self.ptype.wrap_results(r, E)


class Scipy_optimizer_energy(ModuleBase):
    """Optimization module based on scipy.optimize.minimize. Satisfies
    simde.Energy property type. Has two submodules that satisfy
    simde.MolecularBasisSet and simde.AOEnergy property types.
    """
    ptype = Energy()

    def __init__(self):
        ModuleBase.__init__(self)
        self.description(self.__doc__)
        self.satisfies_property_type(self.ptype)
        self.add_submodule(MolecularBasisSet(), "MolecularBasisSet")
        self.add_submodule(AOEnergy(), "AOEnergy")

    def run_(self, inputs, submods):
        system_in, = self.ptype.unwrap_inputs(inputs)
        mol = system_in.molecule
        coords = ct.get_molecule_coordinates(mol)
        symbols = ct.get_molecule_symbols(mol)
        aobasis = submods["MolecularBasisSet"].run_as(MolecularBasisSet(), mol)
        aos = AOSpaceD(aobasis)
        def _get_energy(coords):
            m = ct.get_molecule(symbols, coords)
            chem_sys = ChemicalSystem(m)
            E = submods["AOEnergy"].run_as(AOEnergy(), aos, chem_sys)
            return E
        result = minimize(_get_energy, ravel(coords))
        E = result.fun
        r = self.results()
        return self.ptype.wrap_results(r, E)