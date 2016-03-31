# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:22:00 2016

@author: ko
"""
#midlertiddige klasser
#import classes as cl

#permenente
import numpy as np
import temperature
import salinity
import ice_volume


def singh(world):

    ode = {}    
    ode[world.ocean, 'S'] = world.ocean.salt_func(world)         
    ode[world.ocean, 'T'] = world.ocean.temp_func(world)     
    ode[world.ice, 'V'] = world.ice.V_func(world)
    
    return ode
    
def extra(world, ode):
    return ode
    
    
def singh_gl(world):

    ode = {}
    
    ode[world.ocean, 'S'] = world.ocean.salt_func(world) 
    
    ode[world.ocean, 'T'], ode[world.ice, 'V'] = world.ocean.temp_func(world) 
    return ode
    
    
    
def test(world):
    print 'ode - test'
    ode = {}
    ode[world.ocean, 'test'] = world.ocean.test + 1
    return ode 
    

if __name__ == "__main__":
    import classes as cl    
    import ode_solver as solver   
    
    dim = cl.dim(dt = 1)
    world = cl.world(dim)
 
    world.print_all()   
    world = solver.forward_euler(world)
    world.print_all()    
    
    
#    def forward_test(world,
#                     ode_func = singh,
#                     ode_func_extra = None,
#                     **kw):
#        for i in xrange(world.dim.steps):        
#            ode = ode_func(world)                    
#                    
#            for place in ode:
##                banan = getattr(place[0], place[1]) + ode[place] * world.dim.dt 
#                new = getattr(place[0], place[1]) + ode[place]
#                setattr(place[0], place[1], new)
#            
#        return world
    
    