# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:22:00 2016

@author: Rasmus Nordfang
"""

#permenente
import numpy as np
import temperature
import salinity
import ice_volume



'''
Main ODE Functions. Calculates the ordinary differential functions.
Function pointer to which of the following funtions should be used is stored in world.ode_func


Input: world
Output: Dictionary with key=[placement, paramter-name]  value=dx/dt 
'''

def singh(world): #Default ODE's. 

    ode = {}    
    ode[world.ocean, 'S'] = world.ocean.salt_func(world) #Pointer to salinity function        
    ode[world.ocean, 'T'] = world.ocean.temp_func(world) #pointer to temperature function     
    ode[world.ice, 'V'] = world.ice.V_func(world) #pointer to ice-volume function
    
    return ode
    
def ocean(world): #only ocean
    ode = {}    
    ode[world.ocean, 'S'] = world.ocean.salt_func(world) #Pointer to salinity function         
    ode[world.ocean, 'T'] = world.ocean.temp_func(world) #pointer to temperature function        
    return ode
    
 
 
 
'''
extra ODE functions. Calculates extra ode functions
Function pointer to which funtion should be used is stored in world.ode_func_extra

input: world, dictionary from above functions
output: Dictionary with key=[placement, paramter-name]  value=dx/dt
'''   
    
def extra(world, ode): #Should be empty. This function is used when no extra ODE's are added to the default ones.
    return ode
    



    

if __name__ == "__main__":
    '''
    Only for testing
    '''
    
    import classes as cl    
    import ode_solver as solver   
    
    dim = cl.dim(dt = 1)
    world = cl.world(dim)
 
    world.print_all()   
    world = solver.forward_euler(world)
    world.print_all()    
    
    

    