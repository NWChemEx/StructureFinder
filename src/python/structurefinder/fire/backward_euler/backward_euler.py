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
from simde import (EnergyNuclearGradientStdVectorD, TotalEnergy,TotalEnergyNuclearOptimization, MoleculeFromString)
from berny import Berny, geomlib
import chemist
import numpy as np
import tensorwrapper as tw

class GeomoptViaBackwardEulerFIRE(pp.ModuleBase):

    def __init__(self):
        pp.ModuleBase.__init__(self)
        self.satisfies_property_type(TotalEnergyNuclearOptimization())
        self.description("Performs Backward Euler-FIRE optimization")
        self.add_submodule(TotalEnergy(), "Energy")
        self.add_submodule(EnergyNuclearGradientStdVectorD(), "Gradient")
        self.add_submodule(MoleculeFromString(), "StringConv")

    def run_(self, inputs, submods):
        pt = TotalEnergyNuclearOptimization()
        mol, = pt.unwrap_inputs(inputs)
        molecule = mol.molecule

        # Construction of the coordinate vecotr
        # numb_atoms = molecule.size()  #<-- Number of atoms
        # numb_coord = 3*numb_atoms     #<-- Number of coordinates
        # R_xyz = np.zeros(numb_coord)  #<-- Initializing the coordiantes vector
        # for i in range(molecule.size()):
        #     for j in range(3):
        #         R_xyz[j]   = molecule.at(i).x
        #         R_xyz[j+1] = molecule.at(i).y
        #         R_xyz[j+2] = molecule.at(i).z
        # return R_xyz
        # def e_func(geom):
        #     return submods["Energy"].run_as(TotalEnergy(), geom)
        
        # def grad_func(geom):
        #     return submods["Gradient"].run_as(EnergyNuclearGradientStdVectorD(), geom, geom_points.as_point_set())

        # Loads the geometry string into the Berny optimizer
        # object.
        # optimizer = BE2_FIRE(settings)
        # optimized_energy, optimized_geom = optimizer.optimize(xyz, e_func, grad_func)
        # print(optimized_geom)

       # Optimized energy is of type "float"

#-----------------------------------------------------------------------------------------------------------------
#    def BE2_FIRE(self,v0 = 0,h0 = 0.03, alpha = 0.1, t_max=0.3, numbcycles=1000, error = 10**(-8)):
        
#         #----------------------------------------------------------------------
#         v_initial = np.array(v0)                   #<-- Initial Velocity
#         VEL = [v_initial]                          #<-- Velocity list
#         #----------------------------------------------------------------------
#         h = h0                                     #<-- Time step
#         #----------------------------------------------------------------------
#         #----- Optimization cycle parameters ----------------------------------
#         Ncycle = numbcycles
#         Np = 0
#         Nreset = 0
#         #------FIRE parameters -----------------------------------------------
#         alpha = alpha
#         t_max = t_max
#         error = error
#         #------Counters ------------------------------------------------------
#         k=0                                   #<-- Convergence cycles
#         i=0                                   #<-- numpy array counter
#         #----------------------------------------------------------------------
#         while k < (Ncycle):
#             k = k + 1
#             #------------ NORMALIZED POS ERROR CONVERGENCE --------------------
#             up_pos_norm_error = np.linalg.norm(POS[i] - MIN)
#             down_pos_norm_error = np.linalg.norm(POS[0]-MIN)
#             pos_norm_error = (up_pos_norm_error/down_pos_norm_error)
#             #------------------------------------------------------------------
#             if pos_norm_error < error:
#                 break
#             #----- FIRE ------------------------------------------------------
#             #---- alpha RESET and Half time step setting ---------------------
#             if  np.dot(VEL[i],F[i])<=0:
#                 VEL[i] = v_initial
#                 h = 0.5*h
#                 alpha = 0.1
#                 Np = 0
#                 Nreset = Nreset + 1
                        
#             elif np.dot(VEL[i],F[i])>0:
#                 Np = Np + 1
                
#                 #----- FORCE UNITE VECTOR ------------------------------------
                
#                 fh = F[i]/np.linalg.norm(F[i])
                
#                #---- FIRE VELOCITY CORRECTION -------------------------------
#                 VEL[i] = (1-alpha)*VEL[i] + alpha*np.linalg.norm(VEL[i])*fh
#                #-------------------------------------------------------------
               
#                #----- alpha, time step UPGRADE ------------------------------
#                 if  Np>5:
                                    
#                     h = min(1.1*h,t_max)
                    
                    
#                     alpha = 0.99*alpha
                
                
#             #----- NEW ELEMENT in the lists ----------------------------------- 
#             POS.append(0)
#             F.append(0)
#             norm_force.append(0)
#             VEL.append(0)
#             E.append(0)
#             #--------------- BACKWARD EULER(2)---------------------------------
#             #-- POSITION ------------------------------------------------------
#             POS[i+1] = POS[i] + VEL[i]*h 
#             #-- FORCE ---------------------------------------------------------
#             F[i+1] = Lennard_Jones_Potential(POS[i+1]).force_cart
#             norm_force[i+1] =  np.linalg.norm(F[i+1])
#             #-- VELOCITY ------------------------------------------------------
#             VEL[i+1] = VEL[i] + F[i+1]*h
#             #-- ENERGY -------------------------------------------------------
#             E[i+1] = Lennard_Jones_Potential(POS[i+1]).LJ_energy
                        
#             #-------------- NORMALIZED ENERGY ERROR CONVERGENCE ---------------            
#             if abs(E[i+1]-MIN)/abs(E[0]-MIN) < error:
                
#                 break
#             #------------------------------------------------------------------
            
#             #------------- NORM OF FORCE CONVERGENCE --------------------------
#             if np.linalg.norm(F[i+1])< error:
                
#                 break
#             #------------------------------------------------------------------         

            
#             #--- LIST position update------------------------------------------
#             i = i + 1
#         #--------- END OF OPTIMIZTION PROCESS ---------------------------------
        
#         #--- COLLECTED DATA ANALYSIS ------------------------------------------
#         #----------------------------------------------------------------------
#         POS = np.array(POS) #<--- To calculate the position errors
#         ini_dist = np.linalg.norm(r0 - pos_minima)
#         #---- position error calculation -------------------------------------
#         position_errors = np.linalg.norm((POS - pos_minima), axis =1)/ini_dist 
#         #---------------------------------------------------------------------
       
#         #----------------- Energy Error calculations -------------------------
#         E = np.array(E)#<---- To calculate Energy errors
#         #---------------------------------------------------------------------
#         energy_errors = abs((E - MIN)/(E[0]-MIN))        
#         #------- ATRIBUTES ----------------------------------------------------
#         self.steps = len(E)
#         self.positions = POS
#         self.minima_error = position_errors
#         self.forces = F
#         self.norm_forces = np.array(norm_force)
#         self.velocity = VEL
#         self.energy = E
#         self.energy_errors = energy_errors
#         self.Np = Np
#         self.alpha = alpha
#         self.Nreset = Nreset
#         self.min_force =min(norm_force)
#         #----------------------------------------------------------------------
        

        e = tw.Tensor(np.array(0))
        ps = chemist.PointSetD()
        ps.push_back(chemist.PointD(1.0,2.0,3.0))
        #       print(e)
        rv = self.results()
        return pt.wrap_results(rv, e, ps)

def load_backwardeulerfire_modules(mm):
    mm.add_module("BackwardEulerFire",  GeomoptViaBackwardEulerFIRE())






