# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 14:53:42 2016

@author: Rasmus Nordfang
"""


'''
Ice interactions.

Input: class World
Output: numpy array - Interaction value 

Caution: Do not modify 'world' values in these functions.
'''



#dV/dt
def dV(world):
    dV = (world.atmos.Q_lw - world.atmos.Q_sw - world.ice.x_ice * world.ocean.Q_int[:,:,0]) / (world.ice.p_ice * world.ocean.L_f) + world.ice.dp_ice
    return dV

#Ice grow rate (without I_export)
def singh(world):
    dV = world.ice.dV - world.ice.I_export / (60.0*60.0*24.0*365.0)
    return dV



if __name__ == "__main__":
    '''
    only for testing
    '''    
    
    import classes as cl
        
    dim = cl.dim()
    world = cl.world(dim)