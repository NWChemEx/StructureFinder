#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:30:09 2024

@author: Felix Rojas
"""

# Imported modules
import numpy as np
from LJ_potential import Lennard_Jones_Potential

class NAG_SUT_NDC():
    
    def __init__(self,r0):
        
        #---- DATA INITIALIZAION ---------------------------------------------
        r0 = np.array([r0])
        POS = [r0]                                         #<-- Position list
        VEL = [np.array(0)]
        MIN = -1                                       #<-- Energy minimum
        pos_minima = 2**(1/6)
        #---------------------------------------------------------------------
        E0=Lennard_Jones_Potential(r0).LJ_energy #<-Initial Energy
        E = [E0]
        #---------------------------------------------------------------------
        F0 = Lennard_Jones_Potential(r0).force_cart #<-- Initial Force
        F = [F0]                                           #<-- Forces List
        norm_force = [np.linalg.norm(F[0])]
        #----------------------------------------------------------------------
        Ncycle = 1000
        #---------------------------------------------------------------------
        alpha = 0.001 
        mu = 0.9
        #---------------------------------------------------------------------
        
        k=0
        i=0        
        while k < (Ncycle):
            k = k + 1
            #------------ NORMALIZED POS ERROR CONVERGENCE --------------------
            up_pos_norm_error = np.linalg.norm(POS[i] - MIN)
            down_pos_norm_error = np.linalg.norm(POS[0]-MIN)
            pos_norm_error = (up_pos_norm_error/down_pos_norm_error)
            #------------------------------------------------------------------
            if pos_norm_error < 10**(-8):
                break   
            #----- NEW ELEMENT in the lists ----------------------------------- 
            POS.append(0)
            VEL.append(0)
            F.append(0)
            norm_force.append(0)
            E.append(0)
            #--------------- NAG-SUTSKEVER-- ---------------------------------
            
            #-- VELOCITY ------------------------------------------------------
            VEL[i+1] = mu*VEL[i] + alpha*F[i]
            #-- POSITION ------------------------------------------------------
            POS[i+1] = POS[i] + VEL[i+1]
            #-- FORCE ---------------------------------------------------------
            F[i+1] = Lennard_Jones_Potential(POS[i+1] + mu*VEL[i+1]).force_cart
            norm_force[i+1] =  np.linalg.norm(F[i+1])
            #-- ENERGY -------------------------------------------------------
            E[i+1]=Lennard_Jones_Potential(POS[i+1]).LJ_energy
                        
            #-------------- NORMALIZED ENERGY ERROR CONVERGENCE ---------------            
            if abs(E[i+1]-MIN)/abs(E[0]-MIN) < 10**(-8):
                
                break
            #------------------------------------------------------------------
            
            #------------- NORM OF FORCE CONVERGENCE --------------------------
            if np.linalg.norm(F[i+1])<10**(-8):
                
                break
            #------------------------------------------------------------------
                            
            
            #--- LIST position update------------------------------------------
            i = i + 1
        
        #--------- END OF OPTIMIZTION PROCESS ---------------------------------
        #--- COLLECTED DATA ANALYSIS ------------------------------------------
        #----------------------------------------------------------------------
        POS = np.array(POS) #<--- To calculate the position errors
        ini_dist = np.linalg.norm(r0 - pos_minima)
        #---- position error calculation -------------------------------------
        position_errors = np.linalg.norm((POS - pos_minima), axis =1)/ini_dist 
        #---------------------------------------------------------------------
       
        #----------------- Energy Error calculations -------------------------
        E = np.array(E)#<---- To calculate Energy errors
        #---------------------------------------------------------------------
        energy_errors = abs((E - MIN)/(E[0]-MIN))        
        #------- ATRIBUTES ----------------------------------------------------
        self.steps = len(E)
        self.positions = POS
        self.minima_error = position_errors
        self.forces = F
        self.norm_forces = np.array(norm_force)
        self.velocity = VEL
        self.energy = E
        self.energy_errors = energy_errors
        self.min_force =min(norm_force)
       #----------------------------------------------------------------------
       # END of code ---------------------------------------------------------
