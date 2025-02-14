"""
@author: Felix Rojas
"""

import numpy as np
import plugandplay as pp
from simde import TotalEnergy

class LJ_potential(pp.ModuleBase):
    # Module Construct --------------------------------------------------------
    def __init__(self):
        """
        This module Evaluates the Lennard-Jones 1D potential function (E)
        """
        pp.ModuleBase.__init__(self)
        self.description("Lennard-Jones 1D potential function")
        self.satisfies_property_type(TotalEnergy())
    #--------------------------------------------------------------------------
    
    # Module run_ member function ---------------------------------------------
    def run_(self, inputs):
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
        x0 = pt.unwrap_inputs(inputs)
        #-------------- LENNARD-JONES FUNCTION --------------------------------
        E = lambda x: 4*((1/x**12)-(1/x**6))
        #------------- ANALYTIC FORCE -----------------------------------------        
        DE_x = -24*((2/x0**13)-(1/x0**7))
        FC = -DE_x
        #----------------------------------------------------------------------
        rv = self.results()    
        return pt.wrap_results(rv, E(x0))
    #--------------------------------------------------------------------------
    
def load_Lenard_Jones_potential(mm):
    mm.add_module("Lenard-Jones", LJ_potential())