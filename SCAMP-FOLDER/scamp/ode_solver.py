# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:24:22 2016

@author: Rasmus Nordfang
"""
#import classes
import update
import temperature
import salinity
import copy
from ode import singh

    
'''
Model solvers. The solvers are called with scamp.solver_name(world)

Input: world
Output: world 
'''    
    
def forward_euler(world): #Basic forward euler
    for i in xrange(world.dim.steps):      
        
        world = world.update_func(world, i) #Update all values to current step,        
        ode = world.ode_func(world) #calculate dx/dt for state parameters defined in ode_func
        ode = world.ode_func_extra(world, ode) #calculate extra state paramters defined in ode_func_extra                       

        for place in ode:   #eg place[0]='world.ocean'   place[1] = 'T'   ode[place] = 0.1
            new = getattr(place[0], place[1]) + ode[place] * world.dim.dt_s #T_(i+1) = T + (dT/dt)*dt  
            setattr(place[0], place[1], new) 
        
        if i%world.dim.save_fraction == 0:
#            if i/world.dim.save_fraction%10 == 0:
#                print 'i', i
            
            
            try:                
                world.ocean.p_save[i/world.dim.save_fraction+1,:,:,:] = world.ocean.p[:,:,:] 
                world.ocean.T_save[i/world.dim.save_fraction+1,:,:,:] = world.ocean.T[:,:,:]
                world.ocean.S_save[i/world.dim.save_fraction+1,:,:,:] = world.ocean.S[:,:,:]
                world.ice.V_save[i/world.dim.save_fraction+1,:,:] = world.ice.V[:,:]
                world.ice.h_ice_save[i/world.dim.save_fraction+1,:,:] = world.ice.h_ice[:,:]
            except:
                pass
            
        #Save all values. This could be outcommented if full plots are not desired.
#        world.ocean.p_save[i+1,:,:,:] = world.ocean.p[:,:,:] 
#        world.ocean.T_save[i+1,:,:,:] = world.ocean.T[:,:,:]
#        world.ocean.S_save[i+1,:,:,:] = world.ocean.S[:,:,:]
#        world.ice.V_save[i+1,:,:] = world.ice.V[:,:]
#        world.ice.h_ice_save[i+1,:,:] = world.ice.h_ice[:,:]
                
    return world  
    
 

def runge_kutta_4(world): #Runge kutta 4
    a21 = 1/2.0
    a32 = 1/2.0
    a43 = 1.0
    b1 = 1/6.0
    b2 = 1/3.0
    b3 = 1/3.0
    b4 = 1/6.0
    
    world1 = copy.deepcopy(world)#This is only done once so not time consuming but bad for memory
    world2 = copy.deepcopy(world)
    world3 = copy.deepcopy(world)
    world4 = copy.deepcopy(world)   
 
    
    for i in xrange(world.dim.steps):
        
        world = world.update_func(world,i)
        ode1 = world.ode_func(world, world1) #dT1/dt gemt i ode1
        
        for place in ode1:    
            try:             
                new = getattr(world.ocean, place[1]) + ode1[place] * a21 * world.dim.dt_s
                
            except AttributeError:
                new = getattr(world.ice, place[1]) + ode1[place] * a21 * world.dim.dt_s
            setattr(place[0], place[1], new) #gemmer T + a21 * dT1 i world1
            
            
        world1 = world1.update_func(world1,i)
        ode2 = world1.ode_func(world1, world2) #dT2 gemt i ode2
        
        for place in ode2:    
#            print 'her',  place[0]
            try:             
                new = getattr(world.ocean, place[1]) + ode2[place] * a32 * world.dim.dt_s
            except AttributeError:
                new = getattr(world.ice, place[1]) + ode2[place] * a32 * world.dim.dt_s
            setattr(place[0], place[1], new) #gemmer T + a32 * dT2 i world1
            
            
        world2 = world2.update_func(world2,i)
        ode3 = world1.ode_func(world2, world3) #dT3 gemt i ode3
        
        for place in ode3:    
#            print 'her',  place[0]
            try:             
                new = getattr(world.ocean, place[1]) + ode3[place] * a43 * world.dim.dt_s
            except AttributeError:
                new = getattr(world.ice, place[1]) + ode3[place] * a43 * world.dim.dt_s
            setattr(place[0], place[1], new) #gemmer T + a32 * dT2 i world1
            
        world3 = world3.update_func(world3,i)
        ode4 = world1.ode_func(world3, world4) #dT4 gemt i ode3

        for place in ode1:
            
            try:
                new = getattr(world.ocean, place[1]) + world.dim.dt_s * (b1 * ode1[world1.ocean, place[1]] + 
                                                                         b2 * ode2[world2.ocean, place[1]] + 
                                                                         b3 * ode3[world3.ocean, place[1]] +
                                                                         b4 * ode4[world4.ocean, place[1]] )
                                                        
                setattr(world.ocean, place[1], new)
            except:
                new = getattr(world.ice, place[1]) + world.dim.dt_s * (b1 * ode1[world1.ice, place[1]] + 
                                                                       b2 * ode2[world2.ice, place[1]] + 
                                                                       b3 * ode3[world3.ice, place[1]] +
                                                                       b4 * ode4[world4.ice, place[1]] )
                                                      
                setattr(world.ice, place[1], new)
        
        print world.ocean.T        
        
        world.ocean.T_save[i+1,:,:,:] = world.ocean.T[:,:,:]
        world.ocean.S_save[i+1,:,:,:] = world.ocean.S[:,:,:]
        world.ice.V_save[i+1,:,:] = world.ice.V[:,:]
        
    return world



#for i in range(steps): #Forward euler 
#    print 'season','i', i , '<',0.5 * year / dt
##    RUNGE KUTTA 4
#    a21 = 1/2
#    a32 = 1/2
#    a43 = 1
#    b1 = 1/6
#    b2 = 1/3
#    b3 = 1/3
#    b4 = 1/6
#    #k1 obs diff_func is already multiplied with dt
#    [dT1, dS1, dV1, dV_grow1] = diff_func(T,S,V,pres, h, height, season[i], S_save[0,:,:,0], dV_grow, I_export, dt)
#    #k2
#    [dT2, dS2, dV2, dV_grow2] = diff_func(T+a21*dT1, S+a21*dS1, V+a21*dV1,pres, h, height, season[i], S_save[0,:,:,0], dV_grow, I_export, dt)
#    #k3
#    [dT3, dS3, dV3, dV_grow3] = diff_func(T+a32*dT2, S+a32*dS2, V+a32*dV2,pres, h, height, season[i], S_save[0,:,:,0], dV_grow, I_export, dt)
#    #k4
#    [dT4, dS4, dV4, dV_grow4] = diff_func(T+a43*dT3, S+a43*dS3, V+a43*dV3,pres, h, height, season[i], S_save[0,:,:,0], dV_grow, I_export, dt)
#    
#    T = T + b1 * dT1 + b2 * dT2 + b3 * dT3 + b4 * dT4 
#    S = S + b1 * dS1 + b2 * dS2 + b3 * dS3 + b4 * dS4
#    V = V + b1 * dV1 + b2 * dV2 + b4 * dV3 + b4 * dV4
#    dV_grow = (b1 * dV_grow1 + b2 * dV_grow2 + b3 * dV_grow3 + b4 * dV_grow4) / dt[0] #dV is not one of the ODE variable and is not dt dependent
#
##
##
#def diff_func(T,S,V,pres, h, height, season, S_0, dV, I_export, dt):
#    [NL, NW, layers] = np.shape(T)
#    
#    
#    h_ice = V[:,:] #if V > 0.5 -> h_ice = V and h_ice = 1
#    x_ice = np.ones(shape=(NL, NW)); #percentage ice cover
#    np.putmask(h_ice[:,:], V[:,:] <= 0.5, V[:,:]**0.5 / 2) #if V <= 0.5 then h_ice = sqrt(V) / 2  
#    np.putmask(x_ice[:,:], V[:,:] <= 0.5, 2 * h_ice[:,:]) #if V <= 0.5 then x_ice = 2 * h_ice     
#
#    p = dens(S[:,:,:],T[:,:,:], pres[:])#[kg/m]   (3d [kg/m^3])
#    c_w = heatcap(S[:,:,:],T[:,:,:], pres[:]) #[J/(kg*C)]                            
#    m = height[:] * p[:,:,:]  #[kg]
#    
##    [dS] = salt(S[:,:,:], p[:,:,:], h, x_ice[:,:], I_export, dV, dt, S_0)    
#    [dS] = salt.salt(S[:,:,:], p[:,:,:], h, x_ice[:,:], I_export, dV_grow, S_save[0,:,:,0], W_tot_func = W_tot_func)    
#    
#    
#    
#    [dT, dV, dV_grow] = temp(T[:,:,:], p[:,:,:], h, m[:,:,:], c_w[:,:,:], x_ice[:,:], season, h_ice[:,:], I_export) 
#    
#    return dT*dt[0], dS*dt[0], dV*dt[0], dV_grow*dt[0]


   
#def forward_euler_2(world):
#    for i in xrange(world.dim.steps):        
#        world = world.update_func(world, i)                
#        print 'i', i
#        world.print_all()
#        ode = world.ocean.ode_func(world)                             
##        for place in world.ocean.ode_func(world):
#        for place in ode:    
#            #eg place[0]='world.ocean'   place[1] = 'T'   ode[place] = 0.1
#
#            new = getattr(place[0], place[1]) + ode[place] * world.dim.dt_s #T_(i+1) = T + (dT/dt)*dt  
#            print 'NEW - ', place[1], ' = ', new       
#            setattr(place[0], place[1], new)
#        world.ocean.T = temperature.singh(world)
#        world.ocean.S = salinity.singh(world)
#        world.ice.V = temperature.sing(world)            
#            
#     
#        world.ocean.T_save[i+1,:,:,:] = world.ocean.T[:,:,:]
#        world.ocean.S_save[i+1,:,:,:] = world.ocean.S[:,:,:]
#        world.ice.V_save[i+1,:,:] = world.ice.V[:,:]
#        
#    return world  
    
    
if __name__ == '__main__':
    '''
    Only for testing
    '''
    import update
    import classes as cl
    dim = cl.dim(years = 5)
    
    def test_ode(world, world_save=None):
        #BEGIN should be included in ode files
        if world_save==None:
            world_save = world  
        ode = {}
        #STOP should be included in ode files
        
        ode[world_save.ocean, 'S'] = world.ocean.salt_func(world)         
        ode[world_save.ocean, 'T'] = world.ocean.temp_func(world)     
        ode[world_save.ice, 'V'] = world.ice.V_func(world)
        return ode    
    
    world = cl.world(dim, ode_func = test_ode)
    
    runge_kutta_4(world)
    
    world.plot_all()
    
#    world = cl.world(dim)
#    world = forward_euler(world)
#    world.plot_all()
    
    
#    dim = cl.dim()    
#    world = cl.world(dim, update_func = update.extra)
#    
#    forward_euler(world)