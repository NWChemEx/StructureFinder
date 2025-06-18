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
        mol, points = pt.unwrap_inputs(inputs)
        molecule = mol.molecule

        # Construction of the coordinate vecotr
        numb_atoms = molecule.size()  #<-- Number of atoms
        numb_coord = 3*numb_atoms     #<-- Number of coordinates
        R_xyz = np.zeros(numb_coord)  #<-- Initializing the coordiantes vector
        for i_atom in range(molecule.size()):
            R_xyz[3*i_atom]   = molecule.at(i_atom).x
            R_xyz[3*i_atom + 1] = molecule.at(i_atom).y
            R_xyz[3*i_atom + 2] = molecule.at(i_atom).z
        # print(list(R_xyz))
        # molecule.at(0).x = 110
        # print(molecule)

        def updated_molecule_coord(coordinates,molecule):
            numb_atoms = molecule.size()         
            for i_atom in range(numb_atoms):
                molecule.at(i_atom).x = coordinates[3*i_atom]  
                molecule.at(i_atom).y = coordinates[3*i_atom + 1] 
                molecule.at(i_atom).z = coordinates[3*i_atom + 2]
            return molecule
        
        def e_func(new_coord,molecule):
            current_molecule = updated_molecule_coord(new_coord,molecule)
            return submods["Energy"].run_as(TotalEnergy(), chemist.ChemicalSystem(current_molecule))
    
        def grad_func(new_coord, molecule):
            current_molecule = updated_molecule_coord(new_coord, molecule)
            current_points = current_molecule.nuclei.as_nuclei().charges.point_set.as_point_set()
            return submods["Gradient"].run_as(EnergyNuclearGradientStdVectorD(), chemist.ChemicalSystem(current_molecule), current_points)

        def print_pointset(pointset):
            printout = ' '
            for i in range(pointset.size()):
                printout+= '['
                for j in range(3):
                    printout+= str(round(pointset.at(i).coord(j),10)) + ' '
                printout+= ']'
            print(printout)
        
        # #Loads the geometry string into the Berny optimizer object.
        # optimizer = BE2_FIRE(settings)
        # optimized_energy, optimized_geom = optimizer.optimize(xyz, e_func, grad_func)
        #print(optimized_geom)

    # Optimized energy is of type "float"

#-----------------------------------------------------------------------------------------------------------------
        def backwardeuler_FIRE(molecule,numb_coord, R_xyz, v0 = 0,h0 = 0.03, alpha = 0.1, t_max=0.3, numbcycles=100, error_power = -8, time_step_update = 5):
            mol = molecule
            mol_coord = R_xyz
            #----------------------------------------------------------------------
            error = 10^(error_power)
            v_0 = np.full(numb_coord,v0)         #<-- Initial Velocity
            VEL = [v_0]                          #<-- Velocity vector

            E_0 = e_func(mol_coord,mol)          #<-- Initial Energy
            Energy = [E_0]                       #<-- Energy list
            G_0 =  grad_func(mol_coord,mol)      #<-- Initial Gradient
            Gradient = [G_0]                     #<-- Gradient list
            R_0 = R_xyz                          #<-- Initial Coordinate Vector
            coord_vector = np.array(R_0)         #<-- Coordinate Vector 
            #----------------------------------------------------------------------
            h = h0                                     #<-- Time step
            #----------------------------------------------------------------------
            #----- Optimization cycle parameters ----------------------------------
            Ncycle = numbcycles
            Np = time_step_update
            Nreset = 0
            #------FIRE parameters -----------------------------------------------
            alpha = alpha
            t_max = t_max
            error = error
            #------Counters ------------------------------------------------------
            k=0                                   #<-- Convergence cycles
            i=0                                   #<-- numpy array counter
            #----------------------------------------------------------------------
        backwardeuler_FIRE(molecule,numb_coord, R_xyz)    
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
        

        e = tw.Tensor(np.array([0.0]))
        ps = chemist.PointSetD()
        for i in range(numb_atoms):
            ps.push_back(chemist.PointD(R_xyz[3*i],R_xyz[3*i + 1], R_xyz[3*i + 2]))
        rv = self.results()
        return pt.wrap_results(rv, e, ps)

def load_backwardeulerfire_modules(mm):
    mm.add_module("BackwardEulerFire",  GeomoptViaBackwardEulerFIRE())






