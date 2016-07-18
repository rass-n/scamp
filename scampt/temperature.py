# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 11:04:16 2016

@author: Rasmus Nordfang
"""

#import bohrium as np
import numpy as np
import density
import heat

'''
Temperature interactions.

Input: class World
Output: numpy array - Interaction value 

Caution: Do not modify 'world' values in these functions.
'''


def Q_trans(world): #Temperature from AMOC current
    Q_trans = np.zeros(shape=(world.dim.n_x, world.dim.n_y, world.ocean.layers)) #[W/m^2]                          #Ice layer
    Q_trans[:,:,0] = (1 - world.ice.x_ice) * 2/3.0 * world.ocean.Q_tot_no_ice  #ML
    Q_trans[:,:,1] = (1 - world.ice.x_ice) * 1/3.0 * world.ocean.Q_tot_no_ice  #PC
    Q_trans[:,:,2] = world.ocean.Q_tot_ice * world.ice.x_ice                #DP
    Q_trans[:,:,3] = 0                                #AB    
    return Q_trans
    
def Q_int(world): #Vertical interfacial fluxes
    Q_int = np.zeros(shape=(world.dim.n_x, world.dim.n_y, world.ocean.layers))
    Q_int[:,:,0] = world.ocean.C_0 * (world.ocean.T[:,:,0] - world.ocean.T_bot) #ice - ML
    Q_int[:,:,1] = 2 * world.ocean.k_t_tot[:,:,1] * world.ocean.c_w[:,:,1] * (world.ocean.p[:,:,0] * world.ocean.T[:,:,0] - world.ocean.p[:,:,1] * world.ocean.T[:,:,1]) / (world.ocean.h[0] + world.ocean.h[1]) #ML- PC
    Q_int[:,:,2] = 2 * world.ocean.k_t_tot[:,:,2] * world.ocean.c_w[:,:,2] * (world.ocean.p[:,:,1] * world.ocean.T[:,:,1] - world.ocean.p[:,:,2] * world.ocean.T[:,:,2]) / (world.ocean.h[1] + world.ocean.h[2]) #PC - DP
    Q_int[:,:,3] = 2 * world.ocean.k_t_tot[:,:,3] * world.ocean.c_w[:,:,3] * (world.ocean.p[:,:,2] * world.ocean.T[:,:,2] - world.ocean.p[:,:,3] * world.ocean.T[:,:,3]) / (world.ocean.h[2] + world.ocean.h[3]) #DP - AB
    return Q_int
    
def Q_tot(world): #Total vertical temperature fluxes
    Q_tot=np.zeros(shape=(world.dim.n_x,world.dim.n_y,4)) #[w/m^2]
    Q_tot[:,:,0] = (1 - world.ice.x_ice) * (world.atmos.Q_sw - world.atmos.Q_lw - world.ocean.Q_turb) - world.ice.x_ice * world.ocean.Q_int[:,:,0] - world.ocean.Q_int[:,:,1] + world.ocean.Q_trans[:,:,0]# + Q_int_xy[:,:,0] + Q_int_dia_tot[:,:,0]  #ML  
    Q_tot[:,:,1] = world.ocean.Q_trans[:,:,1] - world.ocean.Q_int[:,:,2] + world.ocean.Q_int[:,:,1]# + Q_int_xy[:,:,1] + Q_int_dia_tot[:,:,1]  #PC
    Q_tot[:,:,2] = world.ocean.Q_trans[:,:,2] - world.ocean.Q_int[:,:,3] + world.ocean.Q_int[:,:,2]# + Q_int_xy[:,:,2] + Q_int_dia_tot[:,:,2] #DP
    Q_tot[:,:,3] = world.ocean.Q_int[:,:,3]# + Q_int_xy[:,:,3] + Q_int_dia_tot[:,:,3]#AB
    return Q_tot
    
        
def Q_turb(world): #Turbulent mixing in top layer
    Q_turb = world.ocean.B_t * (world.ocean.T[:,:,0] - world.ocean.T_ref) #[w/m^2] 
    return Q_turb

def Q_sw(world): #short wave radiation
    #if season==0(winter) -> Q_sw=0,  if season==1(summer) -> Q_sw = the equation
    Q_sw =np.zeros(shape=(world.dim.n_x, world.dim.n_y)) + world.dim.season * (world.ice.x_ice * 80 + (1 - world.ice.x_ice) * 180) #eq 18 Singh
    return Q_sw
   
def Q_lw(world): #long wave radiation
             
             #if season==0         
    Q_lw =((1 - world.dim.season) * ((1 - world.ice.x_ice) * ((world.atmos.wave_A + world.atmos.wave_B * world.ocean.T[:,:,0]) /  world.ocean.n_w - world.atmos.wave_D / 2.0) #over open ocean winter
                                        + world.ice.x_ice  * ((world.atmos.wave_A + world.atmos.wave_B * world.ice.T_top) / world.ocean.n_w - world.atmos.wave_D/2.0))         #over sea ice winter eq 20 Singh
             #if x_ice==0
             + world.dim.season  * ((1 - world.ice.x_ice) * ((world.atmos.wave_A + world.atmos.wave_B * world.ocean.T[:,:,0]) /  world.ocean.n_s - world.atmos.wave_D / 2.0)   #open ocean summer
                                        + world.ice.x_ice  *  (world.atmos.wave_A / world.ocean.n_s - world.atmos.wave_D / 2.0)))                            #sea ice summer eq 20 Singh
    return Q_lw
    
    
def T_top(world): #temperature on the top of the ice-layer
    T_top = np.zeros(shape=(world.dim.n_x, world.dim.n_y))
    T_top[:,:] = (world.ice.x_ice > 0) * (1 / (world.atmos.wave_B/world.ocean.n_w + world.ocean.k_c / world.ice.h_ice[:,:]) * (world.ocean.k_c * world.ocean.T_bot / world.ice.h_ice[:,:] - world.atmos.wave_A / world.ocean.n_w + world.atmos.wave_D / 2.0 ))
    T_top[np.isnan(world.ice.T_top)]=0 #T_top are not be used when T_top=nan but because of problems with mask, it is nescesarry.       
    return T_top
    
def dp_ice_singh(world): #polynya ice added.
    dp_ice = np.zeros(shape=(world.dim.n_x,world.dim.n_y))
    putmask = world.ice.h_ice >= 1
    dp_ice[:,:] = dp_ice + (1 - world.dim.season) * (~putmask * dp_ice + putmask * (2/(60.0*60.0*24.0*365.0))) #[m/s] int_(cold season) dp_ice dt = 1
#    dp_ice[:,:] = dp_ice + (1-world.dim.season) * (~putmask * dp_ice + putmask * (2/(60*60*24*365))) #[m/s] int_(cold season) dp_ice dt = 1    
    return dp_ice
    
def dp_ice(world):
    dp_ice = 0
    return dp_ice    
    
    

def p(world): #Water density each layer
    p = density.dens(world.ocean.S[:,:,:],world.ocean.T[:,:,:], world.ocean.press[:])
    return p

def c_w(world): #Sepcific heat capacity
    c_w = heat.heatcap(world.ocean.S[:,:,:],world.ocean.T[:,:,:], world.ocean.press[:])    
    return c_w

def Q_horiz(world): #Horizontal interfacial fluxes
    
    #x-direction
    Q_int_x = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y, world.ocean.layers)) #[kg * m/s * g/kg]
    Q_int_x[1:world.dim.n_x,:,:] = world.ocean.k_t_horiz_tot[0,0,:] * (world.ocean.p[0:world.dim.n_x-1,:,:] * world.ocean.T[0:world.dim.n_x-1,:,:] - world.ocean.p[1:world.dim.n_x,:,:] * world.ocean.T[1:world.dim.n_x,:,:]) / (world.ocean.h[:])

    #y-direction
    Q_int_y = np.zeros(shape=(world.dim.n_x, world.dim.n_y+1, world.ocean.layers))
    Q_int_y[:,1:world.dim.n_y,:] = world.ocean.k_t_horiz_tot[0,0,:] * (world.ocean.p[:,0:world.dim.n_y-1,:] * world.ocean.T[:,0:world.dim.n_y-1,:] - world.ocean.p[:,1:world.dim.n_y,:] * world.ocean.T[:,1:world.dim.n_y,:]) / (world.ocean.h[:])
    
    #total x and y directions
    Q_int_xy = Q_int_x[0:world.dim.n_x,:,:] - Q_int_x[1:world.dim.n_x+1,:,:] + Q_int_y[:,0:world.dim.n_y,:] - Q_int_y[:,1:world.dim.n_y+1,:]
    
    #diagonal. x=0 y=0 to x=1, y=1
    Q_int_dia_xy = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y+1, world.ocean.layers))
    Q_int_dia_xy[1:world.dim.n_x, 1:world.dim.n_y,:] =  world.ocean.k_t_horiz_tot[0,0,:] * (world.ocean.p[0:world.dim.n_x-1, 0:world.dim.n_y-1,:] * world.ocean.T[0:world.dim.n_x-1, 0:world.dim.n_y-1,:] - world.ocean.p[1:world.dim.n_x, 1:world.dim.n_y,:] * world.ocean.T[1:world.dim.n_x, 1:world.dim.n_y,:] ) / world.ocean.h[:]# * 2**(0.5)

    #diagonal X=0 y=1 to x=1 y=0
    Q_int_dia_yx = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y+1,world.ocean.layers))
    Q_int_dia_yx[1:world.dim.n_x, 1:world.dim.n_y,:] = world.ocean.k_t_horiz_tot[0,0,:] * (world.ocean.p[0:world.dim.n_x-1, 1:world.dim.n_y,:] * world.ocean.T[0:world.dim.n_x-1, 1:world.dim.n_y, :] - world.ocean.p[1:world.dim.n_x, 0:world.dim.n_y-1, :] * world.ocean.T[1:world.dim.n_x, 0:world.dim.n_y-1, :]) / world.ocean.h[:]# * 2 ** (0.5)

    #diagonal total
    Q_int_dia_tot = -Q_int_dia_xy[1:world.dim.n_x+1,1:world.dim.n_y+1,:] -Q_int_dia_yx[1:world.dim.n_x+1,0:world.dim.n_y,:] + Q_int_dia_xy[0:world.dim.n_x, 0:world.dim.n_y,:] + Q_int_dia_yx[0:world.dim.n_x, 1:world.dim.n_y+1,:]

    #total horizontal    
    Q_horiz = Q_int_xy + Q_int_dia_tot / (2.0**(0.5))

    return Q_horiz


def singh(world): #Singh model
    dT = (world.ocean.Q_tot + world.ocean.Q_horiz) / (world.ocean.c_w * world.ocean.m)  #[C/s]    
    
    
    
#    if world.ocean.p[0,0,0] > world.ocean.p[0,0,1]:
#        print 'convection ML - PC'
#        dT[0,0,0] = 0
#        dT[0,0,1] = 0
#        T_tot = (world.ocean.T[0,0,0] * world.ocean.m[0,0,0]  + world.ocean.T[0,0,1] * world.ocean.m[0,0,1]) / (world.ocean.m[0,0,0] + world.ocean.m[0,0,1])
#        world.ocean.T[:,:,0] = T_tot
#        world.ocean.T[:,:,1] = T_tot
#        
#        
#    if world.ocean.p[0,0,1] > world.ocean.p[0,0,2]:
#        print 'Convection PC- DB'
#        dT[0,0,1] = 0
#        dT[0,0,2] = 0
#        T_tot = (world.ocean.T[0,0,1] * world.ocean.m[0,0,1]  + world.ocean.T[0,0,2] * world.ocean.m[0,0,2]) / (world.ocean.m[0,0,1] + world.ocean.m[0,0,2])
#        world.ocean.T[:,:,1] = T_tot
#        world.ocean.T[:,:,2] = T_tot
    
    mask = world.ocean.p[:,:,1] > world.ocean.p[:,:,2]    
    dT[:,:,2] = ~mask*dT[:,:,2] + mask * 0
    dT[:,:,3] = ~mask*dT[:,:,3] + mask * 0    
    T_tot = (world.ocean.T[:,:,1] * world.ocean.m[:,:,1]  + world.ocean.T[:,:,2] * world.ocean.m[:,:,2]) / (world.ocean.m[:,:,1] + world.ocean.m[:,:,2])    
    
    world.ocean.T[:,:,1] = ~mask * world.ocean.T[:,:,1] + mask * T_tot
    world.ocean.T[:,:,2] = ~mask * world.ocean.T[:,:,2] + mask * T_tot   
    
    
#    if world.ocean.p[0,0,2] > world.ocean.p[0,0,3]:
#        print 'Convection DB - AB' 
#        dT[0,0,2] = 0
#        dT[0,0,3] = 0
#        T_tot = (world.ocean.T[0,0,2] * world.ocean.m[0,0,2]  + world.ocean.T[0,0,3] * world.ocean.m[0,0,3]) / (world.ocean.m[0,0,2] + world.ocean.m[0,0,3])
#        world.ocean.T[:,:,2] = T_tot
#        world.ocean.T[:,:,3] = T_tot
    
    return dT 



#old version of Singh.
#def singh_gl(world): 
#    '''
#    INPUT
#    T - temperature [°C]
#    p - density     [kg/m]
#    h - layer depth [m]
#    c_w - specific heat capacity [J / (K * kg)]
#    x_ice - ice parameter. [] number between 0 and 1
#    season [0 winter, 1 summer]
#    h_ice [m]
#    I_export [m] Ice export factor
#
#    OUTPUT
#    T (dT/dt) [C/s]
#    V (dV/dt) [m/s]
#    '''
#
##    [NL, NW, layers] = T.shape
#    #PARAMATERS
#  
#    #OHFC paramaters
#    Q_tot_NoIce = 3 #Unit W/m^2 look at page 4343 middle
#    Q_tot_Ice = 0.5 #Unit W/m^2
#    A = 320 #[W/m^2]
#    B = 4.6 #[w/m^2]
#    D = 90 #[W/m^2]
#    T_ref = -1.8 #[C] minimum temp for ML layer
#    T_bot = -1.8 #[C] ice temp bottom ice present winter eq. 3 Singh
#    B_t = 5.0 #[w/(m^2 * C)]
#    n_w = 2.5
#    n_s = 2.8
#    k_c = 2.0 #[w/(m * C)] thermal conductivity of ice
#    C_0 = 20 #[W/(m^2 * C)]
#    p_ice = 910 #[kg/m] density ice
#    L_f = 297000 #[J/kg] #heat of fusion for seawater. Taken as constant S=4 [0/00] T = -1.8 [C]
#
#    #mixing parameters #Unit m^2/s                          #index      0      1      2      3
#    k_t = np.array([[0, 10**(-4),   10**(-5), 10**(-7)],    #ICE     ice-ML, ML-PC, PC-DP, DP-AB
#                    [0, 6*10**(-4), 10**(-5), 10**(-7)]])   #NO ICE  ice-ML, ML-PC, PC-DP, DP-AB
#
#    #combine both with and without ice
#    
#    k_tot = k_t[None,0,:] * world.ice.x_ice[:,:,None] + (1 - world.ice.x_ice[:,:,None]) * k_t[None,1,:] #OBS dimensionerne kunne være mærkelige efter classes edit
#
#    #mixing paramters horizontal
#    k_horiz = np.array([10**(-4), 10**(-1), 10**(-1), 10**(-1)])
#
#  
#    Q_trans = np.zeros(shape=(world.dim.n_x, world.dim.n_y, world.ocean.layers)) #[W/m^2]                          #Ice layer
#    Q_trans[:,:,0] = (1 - world.ice.x_ice) * 2/3 * Q_tot_NoIce  #ML
#    Q_trans[:,:,1] = (1 - world.ice.x_ice) * 1/3 * Q_tot_NoIce  #PC
#    Q_trans[:,:,2] = Q_tot_Ice * world.ice.x_ice                #DP
#    Q_trans[:,:,3] = 0                                #AB    
#    
#
#    Q_int = np.zeros(shape=(world.dim.n_x, world.dim.n_y, world.ocean.layers))
#    Q_int[:,:,0] = C_0 * (world.ocean.T[:,:,0] - T_bot) #ice - ML
#    Q_int[:,:,1] = 2 * k_tot[:,:,1] * world.ocean.c_w[:,:,1] * (world.ocean.p[:,:,0] * world.ocean.T[:,:,0] - world.ocean.p[:,:,1] * world.ocean.T[:,:,1]) / (world.ocean.h[0] + world.ocean.h[1]) #ML- PC
#    Q_int[:,:,2] = 2 * k_tot[:,:,2] * world.ocean.c_w[:,:,2] * (world.ocean.p[:,:,1] * world.ocean.T[:,:,1] - world.ocean.p[:,:,2] * world.ocean.T[:,:,2]) / (world.ocean.h[1] + world.ocean.h[2]) #PC - DP
#    Q_int[:,:,3] = 2 * k_tot[:,:,3] * world.ocean.c_w[:,:,3] * (world.ocean.p[:,:,2] * world.ocean.T[:,:,2] - world.ocean.p[:,:,3] * world.ocean.T[:,:,3]) / (world.ocean.h[2] + world.ocean.h[3]) #DP - AB
#
#    #Temperature diffusion from neighbors in 'Lenght' direction
##    Q_int_L = np.zeros(shape=(NL+1, NW, layers)) #[kg * m/s * g/kg]
##    Q_int_L[1:NL,:,:] = k_horiz[:] * (p[0:NL-1,:,:] * T[0:NL-1,:,:] - p[1:NL,:,:] * T[1:NL,:,:]) / (h[:])
#    Q_int_x = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y, world.ocean.layers)) #[kg * m/s * g/kg]
#    Q_int_x[1:world.dim.n_x,:,:] = k_horiz[:] * (world.ocean.p[0:world.dim.n_x-1,:,:] * world.ocean.T[0:world.dim.n_x-1,:,:] - world.ocean.p[1:world.dim.n_x,:,:] * world.ocean.T[1:world.dim.n_x,:,:]) / (world.ocean.h[:])
#
#    #Temperature diffusion from neighbors in 'Width' direction
##    Q_int_W = np.zeros(shape=(NL, NW+1, layers)) #[kg * m/s * g/kg]
##    Q_int_W[:,1:NW,:] = k_horiz[:] * (p[:,0:NW-1,:] * T[:,0:NW-1,:] - p[:,1:NW,:] * T[:,1:NW,:]) / (h[:])
#    Q_int_y = np.zeros(shape=(world.dim.n_x, world.dim.n_y+1, world.ocean.layers))
#    Q_int_y[:,1:world.dim.n_y,:] = k_horiz[:] * (world.ocean.p[:,0:world.dim.n_y-1,:] * world.ocean.T[:,0:world.dim.n_y-1,:] - world.ocean.p[:,1:world.dim.n_y,:] * world.ocean.T[:,1:world.dim.n_y,:]) / (world.ocean.h[:])
#    
#
#    #combined diffusion nearest neighbor. Can be put together at W_tot instead
##    Q_int_LW = Q_int_L[0:NL,:,:] - Q_int_L[1:NL+1,:,:] + Q_int_W[:,0:NW,:] - Q_int_W[:,1:NW+1,:]
#    Q_int_xy = Q_int_x[0:world.dim.n_x,:,:] - Q_int_x[1:world.dim.n_x+1,:,:] + Q_int_y[:,0:world.dim.n_y,:] - Q_int_y[:,1:world.dim.n_y+1,:]
#    
#
#    #Temperature diffusion from diagonal neighboors. The mixing coeffecient is probably not correct
##    Q_int_dia_L = np.zeros(shape=(NL+1, NW+1, layers))
##    Q_int_dia_L[1:NL, 1:NW,:] =  k_horiz[:] * (p[0:NL-1, 0:NW-1,:] * T[0:NL-1, 0:NW-1,:] - p[1:NL, 1:NW,:] * T[1:NL, 1:NW,:] ) / h[:]# * 2**(0.5)
#    
#    Q_int_dia_xy = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y+1, world.ocean.layers))
#    Q_int_dia_xy[1:world.dim.n_x, 1:world.dim.n_y,:] =  k_horiz[:] * (world.ocean.p[0:world.dim.n_x-1, 0:world.dim.n_y-1,:] * world.ocean.T[0:world.dim.n_x-1, 0:world.dim.n_y-1,:] - world.ocean.p[1:world.dim.n_x, 1:world.dim.n_y,:] * world.ocean.T[1:world.dim.n_x, 1:world.dim.n_y,:] ) / world.ocean.h[:]# * 2**(0.5)
#
#
##    Q_int_dia_W = np.zeros(shape=(NL+1, NW+1,layers))
##    Q_int_dia_W[1:NL, 1:NW,:] = k_horiz[:] * (p[0:NL-1, 1:NW, :] * T[0:NL-1, 1:NW, :] - p[0:NL-1, 1:NW, :] * T[1:NL, 0:NW-1, :]) / h[:]# * 2 ** (0.5)
#
#    Q_int_dia_yx = np.zeros(shape=(world.dim.n_x+1, world.dim.n_y+1,world.ocean.layers))
#    Q_int_dia_yx[1:world.dim.n_x, 1:world.dim.n_x,:] = k_horiz[:] * (world.ocean.p[0:world.dim.n_x-1, 1:world.dim.n_y,:] * world.ocean.T[0:world.dim.n_x-1, 1:world.dim.n_y, :] - world.ocean.p[0:world.dim.n_x-1, 1:world.dim.n_y, :] * world.ocean.T[1:world.dim.n_x, 0:world.dim.n_y-1, :]) / world.ocean.h[:]# * 2 ** (0.5)
#
#
#    #Combined Temperature diffussion from all diagonal neighbors with a sqrt(2)
##    Q_int_dia_tot = -Q_int_dia_L[1:NL+1,1:NW+1,:] -Q_int_dia_W[1:NL+1,0:NW,:] + Q_int_dia_L[0:NL, 0:NW,:] + Q_int_dia_W[0:NL, 1:NW+1,:]
#    Q_int_dia_tot = -Q_int_dia_xy[1:world.dim.n_x+1,1:world.dim.n_y+1,:] -Q_int_dia_yx[1:world.dim.n_x+1,0:world.dim.n_y,:] + Q_int_dia_xy[0:world.dim.n_x, 0:world.dim.n_y,:] + Q_int_dia_yx[0:world.dim.n_x, 1:world.dim.n_y+1,:]
#
#
#    #Q for the ML layer
#    Q_turb = B_t * (world.ocean.T[:,:,0] - T_ref) #[w/m^2]
#    
#
#    #Q_cond = k_c * (T_top - T_bot) / h_ice
#
#    dp_ice = np.zeros(shape=(world.dim.n_x,world.dim.n_y))
#
##    if season == 0: #winter
##        Q_sw = 0 #[W/m^2]eq 18 Singh
##
##        Q_lw = ((1 - x_ice) * ((A + B * T[:,:,0]) /  n_w - D / 2) #over open ocean winter
##        + x_ice * ((A + B * T_top) / n_w - D/2)) #over sea ice winter eq 20 Singh
##
##        putmask = h_ice >= 1
##        dp_ice[:] = ~putmask * dp_ice + putmask * (2/(60*60*24*365)) #[m/s] int_(cold season) dp_ice dt = 1
##        
##    else: #summer
##        Q_sw = x_ice * 80 + (1 - x_ice) * 180 #eq 18 Singh
##
##        Q_lw = ((1 - x_ice) * ((A + B * T[:,:,0]) /  n_s - D / 2) #open ocean summer
##        + x_ice * (A / n_s - D / 2)) #sea ice summer eq 20 Singh
#
#
#
#
#    if world.dim.season == 0: #winter
#
#        world.ice.T_top[:,:] = (world.ice.x_ice > 0) * (1 / (B/n_w + k_c / world.ice.h_ice[:,:]) * (k_c * T_bot / world.ice.h_ice[:,:] - A / n_w + D / 2 ))        
#        
#        world.ice.T_top[np.isnan(world.ice.T_top)]=0 #T_top are not be used when T_top=nan but because of problems with mask, it is nescesarry.       
#
#        Q_sw = 0 #[W/m^2]eq 18 Singh
#        Q_lw = ((1 - world.ice.x_ice) * ((A + B * world.ocean.T[:,:,0]) /  n_w - D / 2.0) #over open ocean winter
#        + world.ice.x_ice * ((A + B * world.ice.T_top) / n_w - D/2.0)) #over sea ice winter eq 20 Singh
#
#
#        putmask = world.ice.h_ice >= 1
#        dp_ice[:] = ~putmask * dp_ice + putmask * (2/(60*60*24*365)) #[m/s] int_(cold season) dp_ice dt = 1
#        
#
#    else: #summer
#
#        Q_sw = world.ice.x_ice * 80 + (1 - world.ice.x_ice) * 180 #eq 18 Singh
#        Q_lw = ((1 - world.ice.x_ice) * ((A + B * world.ocean.T[:,:,0]) /  n_s - D / 2) #open ocean summer
#        + world.ice.x_ice * (A / n_s - D / 2)) #sea ice summer eq 20 Singh
#
#
#
#    Q_tot=np.zeros(shape=(world.dim.n_x,world.dim.n_y,4)) #[w/m^2]
##    for k in range(2,5):
##        Q_tot[k] = Q_trans[k] + Q_int[k] - Q_int[k-1] #N3 should be moved outside the loop. BUT hoooooow
#
#
##    Q_tot[:,:,0] = (1 - x_ice) * (Q_sw - Q_lw - Q_turb) - x_ice * Q_int[:,:,0] - Q_int[:,:,1] + Q_trans[:,:,0] + Q_int_L[0:NL,:,0] - Q_int_L[1:NL+1,:,0] + Q_int_W[:,0:NW,0] - Q_int_W[:,1:NW+1,0] #ML
##    Q_tot[:,:,1] = Q_trans[:,:,1] - Q_int[:,:,2] + Q_int[:,:,1] + Q_int_L[0:NL,:,1] - Q_int_L[1:NL+1,:,1] + Q_int_W[:,0:NW,1] - Q_int_W[:,1:NW+1,1] #PC
##    Q_tot[:,:,2] = Q_trans[:,:,2] - Q_int[:,:,3] + Q_int[:,:,2] + Q_int_L[0:NL,:,2] - Q_int_L[1:NL+1,:,2] + Q_int_W[:,0:NW,2] - Q_int_W[:,1:NW+1,2]#DP
##    Q_tot[:,:,3] = Q_int[:,:,3] + Q_int_L[0:NL,:,3] - Q_int_L[1:NL+1,:,3] + Q_int_W[:,0:NW,3] - Q_int_W[:,1:NW+1,3]#AB
##
##    Q_tot[:,:,0] = (1 - x_ice) * (Q_sw - Q_lw - Q_turb) - x_ice * Q_int[:,:,0] - Q_int[:,:,1] + Q_trans[:,:,0] + Q_int_LW[:,:,0] + Q_int_dia_tot[:,:,0]  #ML
##    Q_tot[:,:,1] = Q_trans[:,:,1] - Q_int[:,:,2] + Q_int[:,:,1] + Q_int_LW[:,:,1] + Q_int_dia_tot[:,:,1]  #PC
##    Q_tot[:,:,2] = Q_trans[:,:,2] - Q_int[:,:,3] + Q_int[:,:,2] + Q_int_LW[:,:,2] + Q_int_dia_tot[:,:,2] #DP
##    Q_tot[:,:,3] = Q_int[:,:,3] + Q_int_LW[:,:,3] + Q_int_dia_tot[:,:,3]#AB
#    
#
#    Q_tot[:,:,0] = (1 - world.ice.x_ice) * (Q_sw - Q_lw - Q_turb) - world.ice.x_ice * Q_int[:,:,0] - Q_int[:,:,1] + Q_trans[:,:,0] + Q_int_xy[:,:,0] + Q_int_dia_tot[:,:,0]  #ML
#    Q_tot[:,:,1] = Q_trans[:,:,1] - Q_int[:,:,2] + Q_int[:,:,1] + Q_int_xy[:,:,1] + Q_int_dia_tot[:,:,1]  #PC
#    Q_tot[:,:,2] = Q_trans[:,:,2] - Q_int[:,:,3] + Q_int[:,:,2] + Q_int_xy[:,:,2] + Q_int_dia_tot[:,:,2] #DP
#    Q_tot[:,:,3] = Q_int[:,:,3] + Q_int_xy[:,:,3] + Q_int_dia_tot[:,:,3]#AB
#
#
#    dT = Q_tot / (world.ocean.c_w * world.ocean.m)  #[C/s]    
#
#    #this unit does not make sense
##    print 'Først', x_ice * (Q_lw - Q_sw - Q_int[1]),'anden', p_ice * L_f * dp_ice,'Tredje', I_export
#
##    print 'Q_lw - ',Q_lw 
##    print 'Q_sw - ', Q_sw
##    print 'Q_int - ', Q_int
##    print 'p_ice - ', p_ice
##    print 'I_export - ', world.ice.I_export
##    print 'L_f - ', L_f
##    print 'dp_ice - ', dp_ice
##    V = ((world.ice.x_ice * (Q_lw - Q_sw - Q_int[:,:,0]) + p_ice * L_f * dp_ice - world.ice.I_export)) / (p_ice * L_f)# / (60*60*24*365)
##    print 'Q_lw', Q_lw
##    print 'Q_sw', Q_sw
##    print 'dp_ice', dp_ice    
#    
#    
#    world.ice.dV = world.ice.x_ice * (Q_lw - Q_sw - Q_int[:,:,0]) / (p_ice * L_f) + dp_ice #This should be in update func
##    print 'dV', world.ice.dV
##    print 'I_export', world.ice.I_export        
#    
#    
#    V = world.ice.dV - world.ice.I_export / (60*60*24*365)
#
##    V = ((x_ice * (Q_lw - Q_sw - Q_int[1]) - I_export) / (p_ice * L_f) )
#
##    print 'dV/dt', V
#
#    if world.ocean.p[0,0,0] > world.ocean.p[0,0,1]:
##        print 'BANAN'
#        dT[0,0,0] = 0
#        dT[0,0,1] = 0
#        T_tot = (world.ocean.T[0,0,0] * world.ocean.m[0,0,0]  + world.ocean.T[0,0,1] * world.ocean.m[0,0,1]) / (world.ocean.m[0,0,0] + world.ocean.m[0,0,1])
#        world.ocean.T[:,:,0] = T_tot
#        world.ocean.T[:,:,1] = T_tot
#        
#        
#    if world.ocean.p[:,:,1] > world.ocean.p[:,:,2]:
##        print 'APPLE'
#        dT[0,0,1] = 0
#        dT[0,0,2] = 0
#        T_tot = (world.ocean.T[0,0,1] * world.ocean.m[0,0,1]  + world.ocean.T[0,0,2] * world.ocean.m[0,0,2]) / (world.ocean.m[0,0,1] + world.ocean.m[0,0,2])
#        world.ocean.T[:,:,1] = T_tot
#        world.ocean.T[:,:,2] = T_tot
#        
#        
#    if world.ocean.p[:,:,2] > world.ocean.p[:,:,3]:
##        print 'ORANGE' 
#        dT[0,0,2] = 0
#        dT[0,0,3] = 0
#        T_tot = (world.ocean.T[0,0,2] * world.ocean.m[0,0,2]  + world.ocean.T[0,0,3] * world.ocean.m[0,0,3]) / (world.ocean.m[0,0,2] + world.ocean.m[0,0,3])
#        world.ocean.T[:,:,2] = T_tot
#        world.ocean.T[:,:,3] = T_tot
#        
#    return dT, V
    
    
if __name__ == "__main__":
    import classes as cl
    import sys
    sys.path.append('/home/ko/Master/backward/Runge-kutta')
    import matplotlib.pyplot as plt    
    import temp
    
    
    dim = cl.dim()
    dim.season = 0   
    
    world = cl.world(dim)
#    dT, dV = singh(world)
#    lw_summer_ice = np.zeros(118)
#    lw_summer_noice = np.zeros(118)
    lw_winter_ice = np.zeros(100)
#    lw_winter_noice = np.zeros(118)
    T = np.zeros(118)   
    Q_ice = np.zeros(118)
    
    def Q_ml(T):
        Q = 20 * (T + 1.8)
        return Q        
        
#    h_ice = np.zeros(100)
    for i in xrange(-18,100, 1):
        T[i+18] = i*0.1
        Q_ice[i+18] = Q_ml(i*0.1)
#        world.ocean.T[:,:,0] = i*0.1
#        
#        
#        world.ice.x_ice = np.array([[1]])
#        world.dim.season = 1            
#        lw_summer_ice[i+18] = Q_lw(world)
#
#        world.ice.x_ice = np.array([[0]])
#        lw_summer_noice[i+18] = Q_lw(world)
#        
#        world.ice.x_ice = np.array([[1]])
#        world.dim.season = 0
##        lw_winter_ice[i+18] = Q_lw(world)
#        
#        world.ice.x_ice = np.array([[0]])
#        lw_winter_noice[i+18] = Q_lw(world)
#
#
##        world.ice.h_ice = np.array([[i*0.1]])
##        h_ice[i] = i*0.1
##        T[i] = T_top(world)
##    world.ice.h_ice = np.array([[1]])
    plt.plot(T, Q_ice)
#    plt.xlabel('T_ML')
#    plt.plot(T, lw_summer_ice /(world.ice.p_ice * world.ocean.L_f), label='Summer - ice')
#    plt.plot(T, lw_summer_noice /(world.ice.p_ice * world.ocean.L_f), label='Summer - open ocean')
#    plt.plot(T, lw_winter_noice /(world.ice.p_ice * world.ocean.L_f), label='winter - open ocean')
#    
    plt.xlabel('T_ML    [C]')
    plt.ylabel('Q_(ML-ice)   [W m$^{-2}$]')
#    plt.legend(loc='upper left')
    plt.show()
#    
#    for i in xrange(0,100):
#        print i*0.1
#        world.dim.season = 0
#        world.ice.x_ice = np.array([[1]])
#        world.ice.h_ice = np.array([[i*0.1]])
#        world.ice.T_top = T_top(world)
#        h_ice[i] = i*0.1
#        lw_winter_ice[i] = Q_lw(world)
#                
#        
#    plt.plot(h_ice, lw_winter_ice/(world.ocean.L_f * world.ice.p_ice), label='winter - ice')    
#    plt.xlabel('h_ice    [m]')
#    plt.ylabel('Q_lw/(p_ice * L_f)    [m s$^{-1}$]')
#    plt.show()
    
        
    
    
    world.ice.x_ice = np.array([[0]])
    world.dim.season = 1
    world.ocean.T[:,:,0] = -1.8
    
    print 'x_ice', world.ice.x_ice
    print 'season', world.dim.season
    print 'T_ML', world.ocean.T[:,:,0]
    print Q_lw(world)
    
    
#    plt.plot(h_ice, T)
#    plt.xlabel('h_ice [m]')
#    plt.ylabel('T^top [C]')
#    plt.show()
#    print world.ice.h_ice
#    print T_top(world)    
#    
    
#    print 'Ny dT \n', dT
#
##    T, p, h, m, c_w, x_ice,season, h_ice, I_export
#    [dT_gl, a, b] = temp.temp(world.ocean.T, world.ocean.p, world.ocean.h, world.ocean.m, world.ocean.c_w, world.ice.x_ice, world.dim.season, world.ice.h_ice, world.ice.I_export)
#    
#    print 'gl dT \n', dT_gl    
#    
#    
#    print 'Ny dV \n', dV
#    
#    print 'Gl dV \n', a
#    
#    print 'Gl dV 2 \n', b
    
    
    
        