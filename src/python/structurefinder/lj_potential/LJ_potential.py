#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 10:51:27 2023

@author: Felix Rojas
"""

import numpy as np

class Lennard_Jones_Potential():
    """
    This Class Evaluates the Spiral Potential Function, and calculates
    the Force in Polar and Cartesian coordinates. T
    
    The object of this class will have as atributes the force values in each 
    coordinate system and the spiral "energy" value.
    
    Also as a way to check the analytical derivatives, numerical derivatives
    have been calculated.
    """
    
    
    def __init__(self,x_coordinate):
        
        x0 = x_coordinate
        
        
        #-------------- LENNARD-JONES FUNCTION ----------------------------
        E = lambda x: 4*((1/x**12)-(1/x**6))
        
        #-----------------------------------------------------------
        
        #------------- ANALYTIC FORCE ------------------------------
        
        DE_x = -24*((2/x0**13)-(1/x0**7))
        FORCE_CARTESIAN = -DE_x
        
        #self.force_polar = FORCE_POLAR
        self.force_cart  = FORCE_CARTESIAN
        self.LJ_energy = E(x0)
        
        # 3. Numerical FORCE calculation ------------------------------
        #    Using Richardson Extrapolation to get O(h**4) by considering
        #    The central formula approximation to the derivative
        #
        #    f'(x0) = (f(x0+h)-f(x0-h))/(2*h),   O(h**2)
        # --------------------------------------------------------------
        #Ec = lambda x, y: (a-x)**2  + b*(y-x**2)**2
 
        #-----------------------------------------------------------------
        #N1x = lambda h: (Ec(x0+h,y0) - Ec(x0-h,y0))/(2*h)
        #Dx2 = (1/3)*(4*N1x(h/2)- N1x(h))
        #-----------------------------------------------------------------
        #N1y = lambda h: (Ec(x0,y0+h) - Ec(x0,y0-h))/(2*h)
        #Dy2 = (1/3)*(4*N1y(h/2)- N1y(h))
        #----------------------------------------------------------------
        #FORCE_NUMERICAL = -np.array([Dx2,Dy2])
               
        #self.Ec = Ec(x0,y0)
        #self.force_numerical = FORCE_NUMERICAL        
        
        
