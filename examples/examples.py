# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:33:02 2016

@author: Rasmus Nordfang
"""


###########################
# Import relevant modules #
###########################

try: #if the package is installed or placed in same directory
    import scamp as sc #import master-project package

except ImportError: #check upper directory for scamp package.
    try:    
        import sys
        sys.path.append("..")
        import scamp as sc
        
    except: #could not find scamp package
        sys.exit("scamp is not installed or placed in the correct directory")
import numpy as np




###################################################
# Simple example. Default model (stable) #
###################################################

dim = sc.dim(years = 100.0, dt =  1/20.0) #define default dimension of the model except large timesteps and short timeframe
world = sc.world(dim) #define default model (ocean, ice-layer and simple atmosphere)
world = sc.forward_euler(world) #extrapolate the world through time
world.plot_all() #plot results



####################################################################
# Modified model. Original Singh model with polynya ice (unstable) #
####################################################################

#dim_1 = sc.dim(dt = 1/20.0, years = 100.0) #define default dimensions
#world_1 = sc.world(dim_1) #define default world
#
#
#def dp_ice_singh(world): #polynya ice added.    
#    dp_ice = np.zeros(shape=(world.dim.n_x, world.dim.n_y))    
#    putmask = world.ice.h_ice >= 1        
#                 #if winter=>1 summer=>0             If h_ice < 1        if h_ice > = 1        
#    dp_ice[:,:] = dp_ice + (1 - world.dim.season) * (~putmask * dp_ice + putmask * (2 / world.dim.year_s)) #[m/s] int_(cold season) dp_ice dt = 1    
#    return dp_ice
#
#world_1.ice.dp_ice_func = dp_ice_singh #change the default function-pointer to the new function
#
#
#
#world_1 = sc.forward_euler(world_1) #extrapolate the world through time
#world_1.plot_all() #plot results




##################################################
# Model with added parameters and modified ocean #                       
##################################################


#dim_2 = sc.dim(dt = 1/20.0, years = 20.0) #define default dimensions 
#ocean_2 = sc.ocean(dim_2) #define ocean seperate
#ocean_2.fish = ocean_2.new_value(100) #Define new parameter 'fish', and set it to 100 in each grid
#world_2 = sc.world(dim_2, ocean_start = ocean_2) #define world including modified ocean
#
#def fish_func(world, i): #define a function for the new parameter fish
#    world.ocean.fish = world.ocean.fish + 4 #This value is only for the purpose of demonstrating the principle
#    return world
#    
#world_2.update_extra_func = fish_func #include the new parameter-function to the model
#
#
#world_2 = sc.forward_euler(world_2) #extrapolate the world throug time
#
#print world_2.ocean.fish #print the new value
