import pluginplay as pp
import simde
import chemist
import .chemisttools as ct
import numpy as np
from scipy.optimize import minimize

class scipy_optimizer(pp.ModuleBase):
    """Optimization module based on scipy.optimize.minimize"""
    ptype = simde.AOEnergy()

    def __init__(self):
        pp.ModuleBase.__init__(self)
        self.description(self.__doc__)
        self.satisfies_property_type(self.ptype)
        self.add_submodule(self.ptype, "AOEnergy")

    def run_(self, inputs, submods):
        [aos, system_in] = self.ptype.unwrap_inputs(inputs)
        mol = system_in.molecule
        coords = ct.get_coordinates(mol)
        symbols = ct.get_symbols(mol)
        def _get_energy(coords):
            m = ct.get_molecule(symbols, coords)
            chem_sys = chemist.ChemicalSystem(m)
            E = submods["AOEnergy"].run_as(self.ptype, aos, chem_sys)
            return E
        result = minimize(_get_energy, np.ravel(coords))
        E = result.fun
        print('Energy: ', E)
        r = self.results()
        return self.ptype.wrap_results(r, E)