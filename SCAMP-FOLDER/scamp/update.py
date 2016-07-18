# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:06:41 2016

@author: Rasmus Nordfang
"""


import numpy as np
from density import dens
from heat import heatcap

'''
Functions to update all parameters defined by ordinary functions and not by differential equations.

Input: world - class world, i - int, the iteration number from the solver

Output: world
'''

#main update func. Should not be modified.
def main(world, i):
    '''
    Execute update functions for ocean, atmos, ice and extra functions.
    
    'update_func' is a function pointer to the desired update function. 
    Per default the Singh functions are used. Eg. world.ocean.update_func = ocean_singh 
    '''

    world.dim.season = 0 if i*world.dim.dt_y % 1 < world.dim.season_factor else 1    
   
    world = world.ice.update_func(world, i) #ATTENTION. x_ice and h_ice are updatet in this function and should be the first update function
    world = world.ocean.update_func(world, i) #update ocean paramters
    world = world.atmos.update_func(world, i) #update atmosphere paramters
    world = world.update_extra_func(world, i) #update 
    
    return world
    
    
def extra(world, i):
    return world
 

def ice_singh(world,i):
    '''
    Execute all functions related to the ice and update the ice-parameters
    The default functions are defined in temperature.py and salinity.py
    '''
    #WARNING. The order could be very important and variables dependent on ocean or ice variables should be put in dependent. 

    world.ice.V[:,:] = (world.ice.V[:,:] < 0) * 0 + ~(world.ice.V[:,:] < 0) * world.ice.V[:,:] #boundary condition
    
    world.ice.h_ice[:,:] = world.ice.V[:,:] #if V > 0.5 -> h_ice = V and h_ice = 1
    
    putmask = world.ice.V[:,:] <= 0.5
    world.ice.x_ice = np.ones(shape=(world.dim.n_x, world.dim.n_y))

    #if V <= 0.5 then x_ice = 2 * h_ice
    world.ice.h_ice[:,:] = ~putmask * world.ice.h_ice[:,:] + putmask *(world.ice.V[:,:]**0.5 / 2.0)
  
    #if V <= 0.5 then h_ice = sqrt(V) / 2
    world.ice.x_ice[:,:] = ~putmask * world.ice.x_ice[:,:] + putmask *(2 * world.ice.h_ice[:,:])  

    #if dV > 0: export = extra else export = dV * dt_s * C_r + extra
    putmask = world.ice.dV > 0     
    world.ice.export[i%(1.0/world.dim.dt_y),:,:] = ~putmask * 0 + putmask * world.ice.dV * world.dim.dt_s * world.ice.C_r
    world.ice.I_export = np.sum(world.ice.export, axis=0)

    world.ice.dp_ice = world.ice.dp_ice_func(world)
    world.ice.T_top = world.ice.T_top_func(world)    

    return world
    
 
def ocean_singh(world,i):
    '''
    Execute all functions related to the ocean and update the ocean-parameters.
    The default functions are defined in temperature.py and salinity.py
    '''
  
    world.ocean.p[:,:,:] = world.ocean.p_func(world) # density
    world.ocean.c_w[:,:,:] = world.ocean.c_w_func(world) #specific heat capacity of a layer.
    world.ocean.m[:,:,:] = world.ocean.h[:] * world.ocean.p[:,:,:]  #[kg] mass         
        
    #vertical mixing paramter    
    world.ocean.k_s_tot = world.ocean.k_s[None,0,:] * world.ice.x_ice[:,:,None] + (1 - world.ice.x_ice[:,:,None]) * world.ocean.k_s[None,1,:]     
    world.ocean.k_t_tot = world.ocean.k_t[None,0,:] * world.ice.x_ice[:,:,None] + (1 - world.ice.x_ice[:,:,None]) * world.ocean.k_t[None,1,:] 
   
    #horizontal mixing paramter
    world.ocean.k_s_horiz_tot = world.ocean.k_s_horiz[None,0,:] * world.ice.x_ice[:,:,None] + (1 - world.ice.x_ice[:,:,None]) * world.ocean.k_s_horiz[None,1,:]     
    world.ocean.k_t_horiz_tot = world.ocean.k_t_horiz[None,0,:] * world.ice.x_ice[:,:,None] + (1 - world.ice.x_ice[:,:,None]) * world.ocean.k_t_horiz[None,1,:]

    #Salinity related
    world.ocean.W_brine = world.ocean.W_brine_func(world)
    world.ocean.W_melt = world.ocean.W_melt_func(world)
    world.ocean.W_int = world.ocean.W_int_func(world)
    world.ocean.W_export = world.ocean.W_export_func(world)
    world.ocean.W_tot = world.ocean.W_tot_func(world)
    world.ocean.W_horiz = world.ocean.W_horiz_func(world)    
    
    #tempeature related
    world.ocean.Q_trans = world.ocean.Q_trans_func(world)
    world.ocean.Q_int = world.ocean.Q_int_func(world)
    world.ocean.Q_turb = world.ocean.Q_turb_func(world)
    world.ocean.Q_tot = world.ocean.Q_tot_func(world)
    world.ocean.Q_horiz = world.ocean.Q_horiz_func(world)
    
    world.ice.dV = world.ice.dV_func(world)#I should put it another place but it works
        
    return world
    
    
def atmos_singh(world,i):
    '''
    Execute all functions related to the atmosphere and update the atmosphere-parameters
    The default functions are defined in temperature.py and salinity.py
    '''    
    
    world.atmos.Q_sw[:,:] = world.atmos.Q_sw_func(world) #shortwave radiation
    world.atmos.Q_lw[:,:] = world.atmos.Q_lw_func(world) #longwave radiation
    return world
    
def ocean_singh_gl(world,i):

    #only for variables independent of each other
    world.ocean.p[:,:,:] = dens(world.ocean.S[:,:,:],world.ocean.T[:,:,:], world.ocean.press[:])#[kg/m]   (3d [kg/m^3])
    world.ocean.c_w[:,:,:] = heatcap(world.ocean.S[:,:,:],world.ocean.T[:,:,:], world.ocean.press[:]) #[J/(kg*C)]    
#    print 'c_w', world.ocean.c_w                        
    world.ocean.m[:,:,:] = world.ocean.height[:] * world.ocean.p[:,:,:]  #[kg]         
        
    world.ocean.k_s_tot = world.ocean.k_s[None,0,:] * world.ice.x_ice[:,:,None] + (1 - world.ice.x_ice[:,:,None]) * world.ocean.k_s[None,1,:]     

    world.ocean.W_brine = world.ocean.W_brine_func(world)
 
    world.ocean.W_melt = world.ocean.W_melt_func(world)
   
    world.ocean.W_int = world.ocean.W_int_func(world)
  
    world.ocean.W_export = world.ocean.W_export_func(world)

    world.ocean.W_tot = world.ocean.W_tot_func(world)
 
    world.ocean.Q_tot = world.ocean.Q_tot_func(world)
    return world
    


def ice_singh_gl(world,i):
    ''' Old version Ice. Not used in the default model anymore but could be used as a test'''
    #WARNING. The order could be very important and variables dependent on ocean or ice variables should be put in dependent. 

    world.ice.V[:,:] = (world.ice.V[:,:] < 0) * 0 + ~(world.ice.V[:,:] < 0) * world.ice.V[:,:] #boundary condition
    world.ice.h_ice[:,:] = world.ice.V[:,:] #if V > 0.5 -> h_ice = V and h_ice = 1
    
    putmask = world.ice.V[:,:] <= 0.5
    world.ice.x_ice = np.ones(shape=(world.dim.n_x, world.dim.n_y))

    #if V <= 0.5 then x_ice = 2 * h_ice
    world.ice.h_ice[:,:] = ~putmask * world.ice.h_ice[:,:] + putmask *(world.ice.V[:,:]**0.5 / 2.0)
  
    #if V <= 0.5 then h_ice = sqrt(V) / 2
    world.ice.x_ice[:,:] = ~putmask * world.ice.x_ice[:,:] + putmask *(2 * world.ice.h_ice[:,:])  

#    print 'dV', world.ice.dV
#    print 'export', world.ice.export
#    print 'I_export', world.ice.I_export
    putmask = world.ice.dV > 0     
    #if dV > 0: export = extra else export = dV * dt_s * C_r + extra
    world.ice.export[i%(1.0/world.dim.dt_y),:,:] = ~putmask * 0 + putmask * world.ice.dV * world.dim.dt_s * world.ice.C_r
    world.ice.I_export = np.sum(world.ice.export, axis=0)

#    print 'EFTER'
#    print 'dV', world.ice.dV
#    print 'export', world.ice.export
#    print 'I_export', world.ice.I_export
    return world

 
def atmos_singh_gl(world, i):
    return world
    






def singh(world, i):#gammel
    #IMPORTANT-Remember the order so something doesn't affect others with old/new results
    world.dim.season = 0 if i*world.dim.dt % 1 < world.dim.season_factor else 1
                        
    world.ice.h_ice[:,:] = world.ice.V[:,:] #if V > 0.5 -> h_ice = V and h_ice = 1
    putmask = world.ice.V[:,:] <= 0.5
    world.ice.x_ice = np.ones(shape=(world.dim.n_x, world.dim.n_y))
    
    #if V <= 0.5 then x_ice = 2 * h_ice
    world.ice.h_ice[:,:] = ~putmask * world.ice.h_ice[:,:] + putmask *(world.ice.V[:,:]**0.5 / 2.0)
    #if V <= 0.5 then h_ice = sqrt(V) / 2
    world.ice.x_ice[:,:] = ~putmask * world.ice.x_ice[:,:] + putmask *(2 * world.ice.h_ice[:,:])       
    
    world.ocean.p[:,:,:] = dens(world.ocean.S[:,:,:],world.ocean.T[:,:,:], world.ocean.p[:])#[kg/m]   (3d [kg/m^3])
    world.ocean.c_w[:,:,:] = heatcap(world.ocean.S[:,:,:],world.ocean.T[:,:,:], world.ocean.p[:]) #[J/(kg*C)]                            
    world.ocean.m[:,:,:] = world.ocean.height[:] * world.ocean.p[:,:,:]  #[kg]          
    
    return world
#    return 

    
def ocean_singh_meget_gl(world,i):
    #only for variables independent of each other
    world.ocean.p[:,:,:] = dens(world.ocean.S[:,:,:],world.ocean.T[:,:,:], world.ocean.press[:])#[kg/m]   (3d [kg/m^3])
    world.ocean.c_w[:,:,:] = heatcap(world.ocean.S[:,:,:],world.ocean.T[:,:,:], world.ocean.press[:]) #[J/(kg*C)]                            
    world.ocean.m[:,:,:] = world.ocean.height[:] * world.ocean.p[:,:,:]  #[kg]         
        
    return world
    
    
    
    
if __name__ == "__main__":
    import classes as cl
    
    dim = cl.dim()
    world = cl.world(dim)
  
    print 'START'
    world.print_all()
  
    world = main(world, 0)
    
    print "\nSLUT"
    world.print_all()

#
##    world = singh(world, 0)
#    singh(world, 0)
#
#    print "\nSLUT"
#    world.print_all()