# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 14:53:42 2016

@author: ko
"""

def dV(world):
    dV = (world.atmos.Q_lw - world.atmos.Q_sw - world.ice.x_ice * world.ocean.Q_int[:,:,0]) / (world.ice.p_ice * world.ocean.L_f) + world.ice.dp_ice
    return dV


def singh(world):
#    print 'Q_lw', world.atmos.Q_lw
#    print 'Q_sw', world.atmos.Q_sw
#    print 'dp_ice', world.ice.dp_ice
#    print ' '
#    print 'season', world.dim.season
#    print 'naturlig', world.ice.x_ice * (world.atmos.Q_lw - world.atmos.Q_sw - world.ocean.Q_int[:,:,0]) / (world.ice.p_ice * world.ocean.L_f) 
#    print 'polynya ', world.ice.dp_ice
#    print 'I_export', world.ice.I_export / (60.0*60.0*24.0*365.0)
#    world.ice.test = world.ice.test + world.ice.dp_ice * world.dim.dt_s
#    world.ice.dV = world.ice.x_ice * (world.atmos.Q_lw - world.atmos.Q_sw - world.ocean.Q_int[:,:,0]) / (world.ice.p_ice * world.ocean.L_f) + world.ice.dp_ice 

#    world.ice.dV = (world.atmos.Q_lw - world.atmos.Q_sw - world.ice.x_ice * world.ocean.Q_int[:,:,0]) / (world.ice.p_ice * world.ocean.L_f) + world.ice.dp_ice

#    print 'dV', world.ice.dV
#    print 'I_export', world.ice.I_export    
        
    dV = world.ice.dV - world.ice.I_export / (60.0*60.0*24.0*365.0)
    
    return dV


if __name__ == "__main__":
    import classes as cl
        
    dim = cl.dim()
    world = cl.world(dim)