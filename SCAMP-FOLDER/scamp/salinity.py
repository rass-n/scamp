# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 11:40:58 2016

@author: ko
"""

import numpy as np


    
def W_tot(world):

    W_tot = np.zeros(shape=(world.dim.n_x, world.dim.n_y, world.ocean.layers))
    
    W_tot[:,:,0] = -world.ocean.W_int[:,:,1] + world.ocean.W_melt # + W_int_LW[:,:,0] + W_int_dia_tot[:,:,0]
    W_tot[:,:,1] = -world.ocean.W_int[:,:,2] + world.ocean.W_int[:,:,1] + world.ocean.W_brine -0.8 * world.ocean.W_export #+ W_int_LW[:,:,1] + W_int_dia_tot[:,:,1]#PC
    W_tot[:,:,2] = -world.ocean.W_int[:,:,3] + world.ocean.W_int[:,:,2] - 0.2 * world.ocean.W_export# + W_int_LW[:,:,2] + W_int_dia_tot[:,:,2]#DP
    W_tot[:,:,3] = world.ocean.W_int[:,:,3]# + W_int_LW[:,:,3] + W_int_dia_tot[:,:,3] #AB
    return W_tot
        
    

def W_brine(world):
    W_brine = np.zeros(shape=(world.dim.n_x, world.dim.n_y))
    putmask = world.ice.dV[:,:] > 0
    W_brine[:,:] = ~putmask * W_brine[:,:] + putmask *(world.ice.p_ice * world.ocean.S_start * world.ice.dV[:,:] / (1- world.ocean.S_start/1000.0))
    return W_brine
    
def W_melt(world):
    W_melt = np.zeros(shape=(world.dim.n_x, world.dim.n_y))
    putmask = world.ice.dV[:,:] <= 0
    W_melt[:,:] = ~putmask * W_melt[:,:] + putmask *(world.ice.p_ice * world.ocean.S_start * world.ice.dV[:,:] / (1- world.ocean.S_start/1000.0))
    return W_melt
    
    
def W_int(world):
    W_int = np.zeros(shape=(world.dim.n_x,world.dim.n_y, world.ocean.layers)) #[kg * m/s * g/kg]
    W_int[...,1:] = (2 * world.ocean.k_s_tot[...,1:] * (world.ocean.p[...,:world.ocean.layers-1] * world.ocean.S[...,:world.ocean.layers-1] - world.ocean.p[...,1:] * world.ocean.S[...,1:]) 
                    / (world.ocean.h[:world.ocean.layers-1] + world.ocean.h[1:])) #ML- PC, PC-DP, DP-AB
    return W_int
    
def W_export(world):
    W_export = world.ice.p_ice * world.ocean.S_start * world.ice.I_export / (1- world.ocean.S_start/1000) / (60.0*60.0*24.0*365.0) #[kg*m/s g/kg] #pos
    return W_export
    
    
def W_horiz(world):
    W_int_x = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y, world.ocean.layers)) #[kg * m/s * g/kg]
    W_int_x[1:world.dim.n_x,:,:] = world.ocean.k_s_horiz_tot[:,:,:] * (world.ocean.p[0:world.dim.n_x-1,:,:] * world.ocean.S[0:world.dim.n_x-1,:,:] - world.ocean.p[1:world.dim.n_x,:,:] * world.ocean.S[1:world.dim.n_x,:,:]) / (world.ocean.h[:])
    
    W_int_y = np.zeros(shape=(world.dim.n_x, world.dim.n_y+1, world.ocean.layers))
    W_int_y[:,1:world.dim.n_y,:] = world.ocean.k_s.k_horiz_tot[:,:,:] * (world.ocean.p[:,0:world.dim.n_y-1,:] * world.ocean.S[:,0:world.dim.n_y-1,:] - world.ocean.p[:,1:world.dim.n_y,:] * world.ocean.S[:,1:world.dim.n_y,:]) / (world.ocean.h[:])
    
    W_int_xy = W_int_x[0:world.dim.n_x,:,:] - W_int_x[1:world.dim.n_x+1,:,:] + W_int_y[:,0:world.dim.n_y,:] - W_int_y[:,1:world.dim.n_y+1,:]
    
    W_int_dia_xy = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y+1, world.ocean.layers))
    W_int_dia_xy[1:world.dim.n_x, 1:world.dim.n_y,:] =  world.ocean.k_s_horiz_tot[:] * (world.ocean.p[0:world.dim.n_x-1, 0:world.dim.n_y-1,:] * world.ocean.S[0:world.dim.n_x-1, 0:world.dim.n_y-1,:] - world.ocean.p[1:world.dim.n_x, 1:world.dim.n_y,:] * world.ocean.S[1:world.dim.n_x, 1:world.dim.n_y,:] ) / world.ocean.h[:]# * 2**(0.5)

    W_int_dia_yx = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y+1,world.ocean.layers))
    W_int_dia_yx[1:world.dim.n_x, 1:world.dim.n_x,:] = world.ocean.k_s_horiz_tot[:,:,:] * (world.ocean.p[0:world.dim.n_x-1, 1:world.dim.n_y,:] * world.ocean.S[0:world.dim.n_x-1, 1:world.dim.n_y, :] - world.ocean.p[0:world.dim.n_x-1, 1:world.dim.n_y, :] * world.ocean.S[1:world.dim.n_x, 0:world.dim.n_y-1, :]) / world.ocean.h[:]# * 2 ** (0.5)

    W_int_dia_tot = -W_int_dia_xy[1:world.dim.n_x+1,1:world.dim.n_y+1,:] -W_int_dia_yx[1:world.dim.n_x+1,0:world.dim.n_y,:] + W_int_dia_xy[0:world.dim.n_x, 0:world.dim.n_y,:] + W_int_dia_yx[0:world.dim.n_x, 1:world.dim.n_y+1,:]

    W_horiz = W_int_xy + W_int_dia_tot

    return W_horiz
    
    
    
def singh(world):    
    
    dS = (world.ocean.W_tot + world.ocean.W_horiz) / (world.ocean.p * world.ocean.h) #[m/s * g/kg]

    if world.ocean.p[0,0,0] > world.ocean.p[0,0,1]:
#        print 'BANAN'
        dS[0,0,0] = 0
        dS[0,0,1] = 0
        S_tot = (world.ocean.S[0,0,0] * world.ocean.m[0,0,0]  + world.ocean.S[0,0,1] * world.ocean.m[0,0,1]) / (world.ocean.m[0,0,0] + world.ocean.m[0,0,1])
        world.ocean.S[:,:,0] = S_tot
        world.ocean.S[:,:,1] = S_tot
        
        
    if world.ocean.p[0,0,1] > world.ocean.p[0,0,2]:
#        print 'APPLE'
        dS[0,0,1] = 0
        dS[0,0,2] = 0
        S_tot = (world.ocean.S[0,0,1] * world.ocean.m[0,0,1]  + world.ocean.S[0,0,2] * world.ocean.m[0,0,2]) / (world.ocean.m[0,0,1] + world.ocean.m[0,0,2])
        world.ocean.S[:,:,1] = S_tot
        world.ocean.S[:,:,2] = S_tot
        
        
    if world.ocean.p[0,0,2] > world.ocean.p[0,0,3]:
#        print 'ORANGE' 
        dS[0,0,2] = 0
        dS[0,0,3] = 0
        S_tot = (world.ocean.S[0,0,2] * world.ocean.m[0,0,2]  + world.ocean.S[0,0,3] * world.ocean.m[0,0,3]) / (world.ocean.m[0,0,2] + world.ocean.m[0,0,3])
        world.ocean.S[:,:,2] = S_tot
        world.ocean.S[:,:,3] = S_tot
    
    return dS
#def salt(S, 
#         p, 
#         h, 
#         x_ice, 
#         I_export, 
#         dV, 
#         S_start, 
#         W_tot_func=W_tot_func, 
#         **kwargs):
#    '''function that calculates the new salinity level
#    INPUT    
#    S - Salinity [g/kg]
#    p - density  [kg/m]
#    x_ice - percent ice cover []
#    I_export - Yearly Ice removel [m]
#    dV/dt - ice growth [m/s]
#    
#    OUTPUT
#    dS/dt [g/(kg*s)]    
#    '''    
def singh_gl(world):    
#    [NL, NW, layers] = np.shape(S)
#    [NL, NW, layers] = S.shape
    p_ice = 910.0 #[kg/m^3] density ice                       #index     0      1      2      3
    S_start = 34.0
    k_s = np.array([[0, 10**(-4),   10**(-6), 10**(-8)],    #ICE    ice-ML, ML-PC, PC-DP, DP-AB
                    [0, 6*10**(-4), 10**(-6), 10**(-8)]])   #NO ICE ice-ML, ML-PC, PC-DP, DP-AB    
    
    #combine both with and without ice    
    k_tot = k_s[None,0,:] * world.ice.x_ice[:,:,None] + (1 - world.ice.x_ice[:,:,None]) * k_s[None,1,:]    
    

    #Mixing paramters in the horizontal plane NOT CORRECT
#    k_horiz = np.array([10**(-3), 10**(-3), 10**(-3), 10**(-3)])  
    
    W_brine = np.zeros(shape=(world.dim.n_x, world.dim.n_y))
    W_melt = np.zeros(shape=(world.dim.n_x, world.dim.n_y))


#    print dV
#    np.putmask(W_brine[:,:], dV > 0, p_ice * S_start * dV  / (1- S_start/1000))  #included in PC-layer; 34 is initial salinity ML-layer
#    np.putmask(W_melt[:,:], dV <= 0, p_ice * S_start * dV  / (1- S_start/1000)) #included in ML-layer; 34 is initial salinity ML-layer
    
    
    #np.putmask(W_brine[:,:], dV_dt[:,:] > 0, p_ice * 34 * dV_dt[:,:] / (1- 34/1000))
    putmask = world.ice.dV[:,:] > 0

    W_brine[:,:] = ~putmask * W_brine[:,:] + putmask *(p_ice * S_start * world.ice.dV[:,:] / (1- S_start/1000.0))
    #np.putmask(W_melt[:,:], dV_dt[:,:] <= 0, p_ice * 34 * dV_dt[:,:] / (1- 34/1000))
    putmask = world.ice.dV[:,:] <= 0
    W_melt[:,:] = ~putmask * W_melt[:,:] + putmask *(p_ice * S_start * world.ice.dV[:,:] / (1- S_start/1000.0))
    
    
#    W_int = np.zeros(shape=(NL,NW,layers)) #[kg * m/s * g/kg]
#    W_int[:,:,1] = 2 * k_tot[:,:,1] * (p[:,:,0] * S[:,:,0] - p[:,:,1] * S[:,:,1]) / (h[0] + h[1]) #ML- PC    
#    W_int[:,:,2] = 2 * k_tot[:,:,2] * (p[:,:,1] * S[:,:,1] - p[:,:,2] * S[:,:,2]) / (h[1] + h[2]) #PC - DP
#    W_int[:,:,3] = 2 * k_tot[:,:,3] * (p[:,:,2] * S[:,:,2] - p[:,:,3] * S[:,:,3]) / (h[2] + h[3]) #DP - AB

    
    W_int = np.zeros(shape=(world.dim.n_x,world.dim.n_y, world.ocean.layers)) #[kg * m/s * g/kg]
    W_int[...,1:] = (2 * k_tot[...,1:] * (world.ocean.p[...,:world.ocean.layers-1] * world.ocean.S[...,:world.ocean.layers-1] - world.ocean.p[...,1:] * world.ocean.S[...,1:]) 
                    / (world.ocean.h[:world.ocean.layers-1] + world.ocean.h[1:])) #ML- PC, PC-DP, DP-AB

   
    W_export = p_ice * S_start * world.ice.I_export / (1- S_start/1000.0) / (60*60*24*365) #[kg*m/s g/kg] #pos


    W_tot = W_tot_func_gl(world, W_int, W_melt, W_brine, W_export)
#    W_tot = W_tot_func(world)

    dS = W_tot / (world.ocean.p * world.ocean.h) #[m/s * g/kg]

    if world.ocean.p[0,0,0] > world.ocean.p[0,0,1]:
#        print 'BANAN'
        dS[0,0,0] = 0
        dS[0,0,1] = 0
        S_tot = (world.ocean.S[0,0,0] * world.ocean.m[0,0,0]  + world.ocean.S[0,0,1] * world.ocean.m[0,0,1]) / (world.ocean.m[0,0,0] + world.ocean.m[0,0,1])
        world.ocean.S[:,:,0] = S_tot
        world.ocean.S[:,:,1] = S_tot
        
        
    if world.ocean.p[0,0,1] > world.ocean.p[0,0,2]:
#        print 'APPLE'
        dS[0,0,1] = 0
        dS[0,0,2] = 0
        S_tot = (world.ocean.S[0,0,1] * world.ocean.m[0,0,1]  + world.ocean.S[0,0,2] * world.ocean.m[0,0,2]) / (world.ocean.m[0,0,1] + world.ocean.m[0,0,2])
        world.ocean.S[:,:,1] = S_tot
        world.ocean.S[:,:,2] = S_tot
        
        
    if world.ocean.p[0,0,2] > world.ocean.p[0,0,3]:
#        print 'ORANGE' 
        dS[0,0,2] = 0
        dS[0,0,3] = 0
        S_tot = (world.ocean.S[0,0,2] * world.ocean.m[0,0,2]  + world.ocean.S[0,0,3] * world.ocean.m[0,0,3]) / (world.ocean.m[0,0,2] + world.ocean.m[0,0,3])
        world.ocean.S[:,:,2] = S_tot
        world.ocean.S[:,:,3] = S_tot
    
    return dS



#OBS BRUGES MULIGVIS IKKE LÆNGERE   
def salt_W_tot(S, 
               p, 
               h, 
               x_ice, 
               I_export, 
               dV, 
               dt, 
               S_start, 
               W_tot_func = W_tot,
               **kwargs):
    '''function that calculates the new salinity level
    INPUT    
    S - Salinity [g/kg]
    p - density  [kg/m]
    x_ice - percent ice cover []
    I_export - Yearly Ice removel [m]
    dV/dt - ice growth [m/s]
    
    OUTPUT
    dS/dt [g/(kg*s)]    
    '''        
    [NL, NW, layers] = np.shape(S)
    
    
    p_ice = 910 #[kg/m^3] density ice   
    k_s = np.array([[0, 10**(-4),   10**(-6), 10**(-8)],    #ice ML-PC, PC-DP, DP-AB
                    [0, 6*10**(-4), 10**(-6), 10**(-8)]])   #no ice ML-PC, PC-DP, DP-AB    
    
    #combine both with and without ice    
    k_tot = k_s[None,0,:] * x_ice[:,:,None] + (1 - x_ice[:,:,None]) * k_s[None,1,:]    
    

    #Mixing paramters in the horizontal plane NOT CORRECT
#    k_horiz = np.array([10**(-3), 10**(-3), 10**(-3), 10**(-3)])  
    
    
    W_brine = np.zeros(shape=(NL, NW))
    W_melt = np.zeros(shape=(NL, NW))

    np.putmask(W_brine[:,:], dV > 0, p_ice * S_start * dV  / (1- S_start/1000))  #included in PC-layer; 34 is initial salinity ML-layer
    np.putmask(W_melt[:,:], dV <= 0, p_ice * S_start * dV  / (1- S_start/1000)) #included in ML-layer; 34 is initial salinity ML-layer
    
    
#    W_int = np.zeros(shape=(NL,NW,layers)) #[kg * m/s * g/kg]
#    W_int[:,:,1] = 2 * k_tot[:,:,1] * (p[:,:,0] * S[:,:,0] - p[:,:,1] * S[:,:,1]) / (h[0] + h[1]) #ML- PC    
#    W_int[:,:,2] = 2 * k_tot[:,:,2] * (p[:,:,1] * S[:,:,1] - p[:,:,2] * S[:,:,2]) / (h[1] + h[2]) #PC - DP
#    W_int[:,:,3] = 2 * k_tot[:,:,3] * (p[:,:,2] * S[:,:,2] - p[:,:,3] * S[:,:,3]) / (h[2] + h[3]) #DP - AB

    
    W_int = np.zeros(shape=(NL,NW, layers)) #[kg * m/s * g/kg]
    W_int[...,1:] = (2 * k_tot[...,1:] * (p[...,:layers-1] * S[...,:layers-1] - p[...,1:] * S[...,1:]) 
                    / (h[:layers-1] + h[1:])) #ML- PC, PC-DP, DP-AB

   
#    W_int = int_flux.W_int(NL, NW, layers, )   
   
   
   
    W_export = p_ice * S_start * I_export / (1- S_start/1000) / (60*60*24*365) #[kg*m/s g/kg] #pos
    
    
    W_tot = np.zeros(shape=(NL, NW, layers))
   
    
    W_tot[:,:,0] =  -W_int[:,:,1] + W_melt # + W_int_LW[:,:,0] + W_int_dia_tot[:,:,0]
    W_tot[:,:,1] =  -W_int[:,:,2] + W_int[:,:,1] + W_brine -0.8 * W_export #+ W_int_LW[:,:,1] + W_int_dia_tot[:,:,1]#PC
    W_tot[:,:,2] =  -W_int[:,:,3] + W_int[:,:,2] -  0.2*W_export# + W_int_LW[:,:,2] + W_int_dia_tot[:,:,2]#DP
    W_tot[:,:,3] =  W_int[:,:,3]# + W_int_LW[:,:,3] + W_int_dia_tot[:,:,3] #AB
        

    W_tot = None

    S = W_tot / (p * h) #[m/s * g/kg]
    
    return S
    
if __name__ == '__main__':
    import sys
    sys.path.append('/home/ko/Master/backward/Runge-kutta')    
    import salt    
    
    import classes as cl
    dim = cl.dim()
    world = cl.world(dim)
    dS = singh(world)
    
    print 'Ny dS \n', dS
    
    dS_gl = salt.salt(world.ocean.S, world.ocean.p, world.ocean.h, world.ice.x_ice, world.ice.I_export, world.ice.dV, 34)
    
    print 'Gamml dS \n', dS_gl