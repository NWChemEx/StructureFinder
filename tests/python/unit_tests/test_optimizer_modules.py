import structurefinder.chemist_tools as ct
from structurefinder.optimizer_modules import scipy_optimizer
import unittest
import nwchemex
import pluginplay as pp
import simde
import chemist

class Test_scipy_optimizer(unittest.TestCase):
    def setUp(self):
        self.mol = ct.get_molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
        self.mm = pp.ModuleManager()
        nwchemex.load_modules(self.mm)
        mod_key = 'scipy_optimizer'
        self.mm.add_module(mod_key, scipy_optimizer())

    def test_scf(self):
        self.mm.change_submod('scipy_optimizer', 'AOEnergy', 'SCF Energy')
        bs = self.mm.at("sto-3g").run_as(simde.MolecularBasisSet(), self.mol)
        aos = chemist.AOSpaceD(bs)
        cs1 = chemist.ChemicalSystem(self.mol)
        energy = self.mm.at("scipy_optimizer").run_as(simde.AOEnergy(), aos, cs1)
        self.assertAlmostEqual(energy, -1.0787343257869195, places=5)
