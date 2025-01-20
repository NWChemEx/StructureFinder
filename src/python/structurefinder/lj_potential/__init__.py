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
import numpy as np
import plugandplay as pp
from simde import TotalEnergy

class LJ_potential(pp.ModuleBase):
    # Module Construct --------------------------------------------------------
    def __init__(self):
        """
        This module Evaluates the Lennard-Jones 1D potential function (E), and 
        calculates minus the gradient (Force) in Cartesian coordinates, 
        according to the relation
        
                                 FC = - dE/dx       
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
        inputs : x-coordinate, 
        TYPE ---> Float
    
        Returns
        -------
        E: Lennard-Jonnes 1D potential Energy, 
        TYPE ---> Float 
        
        FC: Force in cartesian coordinates evaluated at the given input, 
        acording to the relation
                                 FC = - dE/dx  
        
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
        return pt.wrap_results(rv, E,FC)
    #--------------------------------------------------------------------------
    
def load_Lenard_Jones_potential(mm):
    mm.add_module("Lenard-Jones", LJ_potential())
