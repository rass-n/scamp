# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 22:26:40 2016

@author: Rasmus Nordfang
"""

#SCAMP!!!!
#from ode import singh
#from ode_solver import forward_euler
import numpy as np
#from ode import singh
import ode
import temperature
import salinity
import ice_volume
import update
import warnings

import density
import heat
import test
from ani import ani


import matplotlib.pyplot as plt

#import ode_solver
#import bohrium as np

#from ode import singh


class dim(object):
    '''Dimensions
    
    Attributes
    ----------        
    n_x : int 
        Number of grids, x-direction
     
    n_y : int
        Number of grids, y-direction
 
    l_x : float
        Length of grid, x-direction
 
    l_y : float
        Length of grid, y-direction
     
    years : int
        Number of years
    
    dt_y : float
        timestep defined in years
    
    dt_s : float
        timestep in seconds
    
    steps : int
        Number of timesteps
        
    season_factor : float
        Fraction of a year it is winter. Float between 0 - 1 
    
    season : bool
        0 = winter, 1 = Summer
    
    '''
    def __init__(self,
                 n_x = 1, #number of grids x-direction
                 n_y = 1, #number of grids y - direction
                 l_x = 1, #length of a grid x-direction
                 l_y = 1, #length of grid y-direction
                 years = 1.0, #number of years
                 dt = 1.0 / 20.0, #dt in years
                 season_factor = 0.5, #fraction of a year with winter
                 save_fraction = 1, 
                 **kw):
        self.n_x = n_x
        self.n_y = n_y
        self.years = years
        self.year_s = 365.0 * 24.0 * 60.0 * 60.0
        self.dt_y = dt 
        self.dt_s = dt * self.year_s
        self.steps = int(years / float(dt))
        self.season_factor = season_factor     
        self.season = 0 #0 winter, 1 summer
        self.save_fraction = save_fraction
        self.save_steps = int(self.steps / float(self.save_fraction))
               
            
class world(object):
    '''
    Attributes
    ----------   
    dim : class dim
        Dimensions
    
    atmos : class atmos
        Atmosphere    
    
    ice : class ice
        Ice-layer    
    
    ocean : class ocean     
        Ocean
        
    Function pointer
    ----------------
    update_func : <type 'function'>
        Pointer to a function that updates all variables. 
        
    ode_func : <type 'function'>
        Pointer to a function containing all differential functions.
        
    Class functions
    ---------------
    
    '''
    def __init__(self,
                 dim, #class dim
                 ocean_start = None, #class ocean
                 ice_start = None,  #class ice
                 atmos_start = None, #class atmos
                 update_func = update.main,
                 update_extra_func = update.extra,
                 ode_func = ode.singh,
                 ode_func_extra = ode.extra,
                 **kw):

        self.dim = dim
        
        if ocean_start == None:
            self.ocean = ocean(dim)
        else:
            self.ocean = ocean_start
            
        if ice_start == None:
            self.ice = ice(dim)
        else:
            self.ice = ice_start
            
        if atmos_start == None:
            self.atmos = atmos(dim)
        else:
            self.atmos = atmos_start
            
        #FUNCTIONS
        self.update_func = update_func #this is not ideal to change but it is possible if it is desired not to update part of the model.
        
        self.update_extra_func = update_extra_func        
        
        self.ode_func = ode_func
        self.ode_func_extra = ode_func_extra
        
    def print_all(self,
                  start=0,
                  slut=1,
                  **kw):
#        print 'P', self.ocean.p[start:slut]        
        print 'T \n', self.ocean.T, '\n'
        print 'S \n', self.ocean.S, '\n'
        print 'V \n', self.ice.V, '\n'
        return
        
    def save(self, name):
        np.savez(name, T = self.ocean.T_save, S = self.ocean.S_save, V = self.ice.V_save, p = self.ocean.p_save)
        
    def load(self, name):
        npzfile = np.load(name+'.npz')
        self.ocean.T_save = npzfile['T']
        self.ocean.S_save = npzfile['S']
        self.ocean.V_save = npzfile['V']
        self.ocean.p_save = npzfile['p']
        
    
    def print_layer(self, layer):
        if not type(layer) == int:
            if layer == 'ML':
                layer = 0
            elif layer == 'PC':
                layer = 1
            elif layer == 'DP':
                layer = 2
            elif layer == 'AB':
                layer = 3
            elif layer == 'ice':
                print 'V \n', self.ice.V[:,:], '\n'
                return
            else:
                warnings.warn('layer is not correct')
                return
        print 'T \n', self.ocean.T[:,:,layer], '\n'
        print 'S \n', self.ocean.S[:,:,layer], '\n'
        return
     
    def animation(self,
                  para,
                  layer=0,
                  vmin = 30,
                  vmax = 35,                
                  **kw):
        [steps,_,_] = para.shape               
        ani(para[:,:,:], steps-1, vmin = vmin, vmax=vmax)
        return
            
    def print_functions(self):
        '''print an overview of world'''
        print 'This is not done yet, but is going to be awesome'
        
    def plot_all(self):
        
        time = np.arange(len(self.ocean.T_save[:,0,0,0])) * self.dim.dt_y * self.dim.save_fraction
        
        plt.plot(time,self.ocean.S_save[:,0,0,0], 'r', label='ML'), plt.plot(time,self.ocean.S_save[:,0,0,1],'b', label='PC') , plt.plot(time, self.ocean.S_save[:,0,0,2],'c', label='DP')
        plt.plot(time,self.ocean.S_save[:,0,0,3],'g', label='AB')
        plt.xlabel('time (years)')
        plt.ylabel('Salinity [g/kg]' )
        plt.title('Salinity Layers')
        plt.grid(True)
        plt.legend(loc='upper left', shadow=True, fontsize='x-large')
        plt.show()
        
        plt.plot(time,self.ocean.T_save[:,0,0,0], 'r', label='ML'), plt.plot(time,self.ocean.T_save[:,0,0,1],'b', label='PC') , plt.plot(time,self.ocean.T_save[:,0,0,2],'c', label='DP')
        plt.plot(time,self.ocean.T_save[:,0,0,3],'g', label='AB')
        plt.xlabel('time (years)')
        plt.ylabel('Temperature [C]')
        plt.title('Temperature')
        plt.grid(True)
        plt.legend(loc='upper left', shadow=True, fontsize='x-large')
        plt.show()        
            
#        plt.plot(time,self.ocean.p_save[:,0,0,0], 'r', label='ML'), plt.plot(time,self.ocean.p_save[:,0,0,1],'b', label='PC') , plt.plot(time,self.ocean.p_save[:,0,0,2],'c', label='DP')
#        plt.plot(time,self.ocean.p_save[:,0,0,3],'g', label='AB')
#        plt.xlabel('time (years)')
#        plt.ylabel('density []')
#        plt.title('density')
#        plt.grid(True)
#        plt.legend(loc='upper left', shadow=True, fontsize='x-large')
#        plt.show()                   
                        
            
#        plt.plot(time,self.ice.V_save[:,0,0], 'r', label='ML'), plt.plot(time,self.ice.V_save[:,0,0],'b', label='PC') , plt.plot(time,self.ice.V_save[:,0,0],'c', label='DP')
#        plt.plot(time,self.ice.V_save[:,0,0],'g', label='AB')
#        plt.xlabel('time (years)')
#        plt.ylabel('ice_vol [m^3]')
#        plt.title('ice_volume')
#        plt.grid(True)
#        plt.legend(loc='upper left', shadow=True, fontsize='x-large')
#        plt.show()
#        
        plt.plot(time,self.ice.h_ice_save[:,0,0], 'r', label='ML'), plt.plot(time,self.ice.h_ice_save[:,0,0],'b', label='PC') , plt.plot(time,self.ice.h_ice_save[:,0,0],'c', label='DP')
        plt.plot(time,self.ice.h_ice_save[:,0,0],'g', label='AB')
        plt.xlabel('time (years)')
        plt.ylabel('ice height [m]')
        plt.title('ice_height')
        plt.grid(True)
        plt.legend(loc='upper left', shadow=True, fontsize='x-large')
        plt.show()
        
        return        
        
    def plot_con(self):
        sum_year = np.zeros(self.dim.years)
        steps_year = int(1/self.dim.dt_y)
        for i in xrange(self.dim.years):
            for k in xrange(steps_year):
#        print 20.0*i + k
                sum_year[i] += sum(self.ocean.S_save[steps_year*i+k+1,0,0,:] * self.ocean.p_save[steps_year*i+k+1,0,0,:] * self.ocean.h)
        
#        print sum_year
        
        print sum_year[2] / sum_year[-1] 
        plt.plot(sum_year)
        plt.show()
        return
    
    def plot_con_2(self):
        time = np.arange(len(self.ocean.S_save[:,0,0,0])-1) * self.dim.dt_y
        
        sigma = np.zeros(self.dim.steps)
        print sigma.shape
        print time.shape
        for i in xrange(self.dim.steps):
            sigma[i] = sum(self.ocean.S_save[i+1,0,0,:] * self.ocean.p_save[i+1,0,0,:] * self.ocean.h)
                      
#        print sigma[0:10] /  (10.0**8)
#        print sigma[-1:]           
        
#        sigma2 = sigma / (10.0**8)
        
#        sigma2[5] = 1.5
        plt.plot(time, sigma)
        plt.plot(time, np.ones(self.dim.steps) * sigma[2])
        plt.show()
        return
                  
class atmos(object):
    '''atmos
    
    Attributes
    ----------
    
    dim : class dim
        Option to specify specific         
        
    Q_lw : numpy array, shape=(dim.n_x, dim.n_y)
        Longwave heat flux
   
    Q_sw : numy array, shape=(dim.n_x, dim.n_y)
        Shortwave heat flux    
    
    wave_A : float 
        OLR constant. Used in Q_lW_func
        
    wave_B : float
        Temperature ORL constant. Used in Q_lw_func 
    
    wave_D : float
        Atmospheric heat transport flux convergence
        
    Functions
    ---------
    update_func : <type 'function'>
        Collection of all ordinary atmosphere functions that should be updated in each timestep.     
        
    Q_sw_func : <type 'function'>
        Shortwave function

    
    Q_lw_func : <type 'function'>
        Longwave function
    
    '''
    def __init__(self,
                 dim,
                 Q_lw_start  = 42.0, #this number not used
                 Q_sw_start =42.0, #this number not used
                 wave_A = 320.0,
                 wave_B = 4.6,
                 wave_D = 90.0,
                 update_func = update.atmos_singh,
                 Q_sw_func = temperature.Q_sw,
                 Q_lw_func = temperature.Q_lw,
                 **kw):
          
        self.dim = dim
        self.Q_lw =  np.ones(shape=(dim.n_x, dim.n_y)) * Q_lw_start
        self.Q_sw =  np.ones(shape=(dim.n_x, dim.n_y)) * Q_sw_start
        self.update_func = update_func 
        self.Q_sw_func = Q_sw_func
        self.Q_lw_func = Q_lw_func
        self.wave_A = wave_A
        self.wave_B = wave_B
        self.wave_D = wave_D

class ocean(object):
    '''Ocean
    
    Attributes
    ----------
    layers : int
        Number of ocean layers
        
    h : numpy array, shape=(layers)
        Depth of each ocean-layer [m].
        
    press : numpy array, shape=(layers)
        Average pressure at the depth [kg/m^3]
    
    dim : class dim
        If special dimensions are needed for Ocean. Use with caution. 
        
    p : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Water density
        
    c_w : numpy array, shape=(dim.n_x, dim.n_y, layers) 
        Heat capacity of ocean water.
        
    m : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Mass grid [kg]
            
    T : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Ocean Temperature [C]    
    
    S : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Ocean Salinity [g/kg]    
    
    p_save : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Density     
    
    T_save : numpy array, shape=(steps, dim.n_x, dim.n_y, layers)
        Saved array            
    
    S_save : numpy array, shape=(steps, dim.n_x, dim.n_y, layers)
        Saved array
        
    k_s : numpy array, shape=(dim.n_x, dim.n_y, layers)
        salinity flow vertical
    
    k_s_horiz : numpy array, shape=(dim.n_x, dim.n_y, layers)
        salintiy flow horizontal
                          
    Q_tot_no_ice : float
        Total heat flux transport  with no ice  
    
    Q_tot_ice : floar
        Total heat flux transport with ice    
    
    T_ref : float
        A reference temperature set as constant.
        
    T_bot : float
        Temperature bottom of ice
        k_s_tot : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Combined mixing with/wihout salinity temperature       
        
    k_t_tot : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Combined mixing with/wihout ice temperature         
         
    k_s_horiz_tot : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Salinity mixing paramter horizontal    
    
    k_t_horiz_tot : numpy array, shape=(dim.n_x, dim.n_y, layers)
        Temperature mixing parameter horizontal        
        
    W_brine : numpy array, shape=(dim.n_x, dim.n_y)
        Brine salinity flux   
    
    W_melt : numpy array, shape=(dim.n_x, dim.n_y)
        1    
    
    W_int : numpy array, shape=(dim.n_x, dim.n_y, layers)  
        1    
    
    W_horiz : numpy array, shape=(dim.n_x, dim.n_y, layers)
        1            
    
    Q_tot : numpy array, shape=(dim.n_x, dim.n_y, layers)
        1    
    
    Q_trans : numpy array, shape=(dim.n_x, dim.n_y, layers)
        #[W/m^2]        
    
    Q_int : numpy array, shape=(dim.n_x, dim.n_y, layers)
        1    
    
    Q_turb : numpy array, shape=(dim.n_x, dim.n_y, layers)
        1    
    
    Q_horiz : numpy array, shape=(dim.n_x, dim.n_y, layers)
        1
        
    B_t : float
        1
        
    n_w : float
        1        
        
    n_s : float
        1    
    
    k_c : float
        1        
        
    C_0 : float
        1    
    
    k_t : numpy array, shape=(2, layers)
        1    
    
    L_f : float
        Heat of 
        
    Functions
    ---------
    
    W_export_func : update
        apple
    
    '''
    def __init__(self, 
                 dim = dim(),
                 layers = 4,        #ML    PC    DP    AB 
                 T_start = np.array([-1.0  , -1.0  ,  4.0, 1.0]),
                 S_start = np.array([32.0, 32.0, 32.0, 35.0]),
                 h = np.array([50.0, 350.0, 800.0, 3000.0]), #technical depth but 'h' is used in Singh
                 temp_func = temperature.singh,
                 salt_func = salinity.singh,
                 update_func = update.ocean_singh,
                 W_brine_func = salinity.W_brine,
                 W_melt_func = salinity.W_melt,
                 W_int_func = salinity.W_int,
                 W_export_func = salinity.W_export,
                 W_tot_func = salinity.W_tot,
                 W_horiz_func = salinity.W_horiz,
                 Q_trans_func = temperature.Q_trans,
                 Q_int_func = temperature.Q_int,
                 Q_turb_func = temperature.Q_turb,
                 Q_tot_func = temperature.Q_tot,
                 Q_horiz_func = temperature.Q_horiz,
                 p_func = temperature.p,
                 c_w_func = temperature.c_w,
                 **kw):
                          
        self.layers = layers
        self.h = h 
#        self.height = np.array([50.0, 300.0, 400.0, 1800.0]) 
        self.press = np.array([60.4, 360.0, 814.0, 3030.0]) #[kg/m^3] general pressure
        self.dim = dim
        self.p = np.ones(shape=(dim.n_x, dim.n_y, self.layers)) * 1 #NOT CORRECT updated before update.py
        self.c_w = np.ones(shape=(dim.n_x, dim.n_y, self.layers)) * 4000 #NOT CORRECT 
        self.m = np.ones(shape=(dim.n_x, dim.n_y, self.layers)) * 1 #NOT CORRECT     
        
        self.T = np.ones(shape=(dim.n_x, dim.n_y, layers)) * T_start
        self.S = np.ones(shape=(dim.n_x, dim.n_y, layers)) * S_start

        self.p_save = np.zeros(shape=(dim.save_steps+1, dim.n_x, dim.n_y, layers)); self.p_save[0,:,:,:] = self.p
        self.T_save = np.zeros(shape=(dim.save_steps+1, dim.n_x, dim.n_y, layers)); self.T_save[0,:,:,:] = self.T
        self.S_save = np.zeros(shape=(dim.save_steps+1, dim.n_x, dim.n_y, layers)); self.S_save[0,:,:,:] = self.S
        self.S_start = S_start[0] #ML start temp. Used in ice
        
        #Constants
        self.k_s = np.array([[0.0, 10.0**(-4),   10.0**(-6), 10.0**(-8)],    #ICE    ice-ML, ML-PC, PC-DP, DP-AB
                             [0.0, 6.0*10**(-4), 10.0**(-6), 10.0**(-8)]])   #NO ICE ice-ML, ML-PC, PC-DP, DP-AB  
                             
        self.k_s_horiz =1000 * np.array([[10.0**(-4), 10.0**(-4),   10.0**(-6), 10.0**(-8)],    #ICE    ice-ML, ML-PC, PC-DP, DP-AB
                                   [6.0*10**(-4), 6.0*10**(-4), 10.0**(-6), 10.0**(-8)]])   #NO ICE ice-ML, ML-PC, PC-DP, DP-AB 
                                   
        
                        
        self.Q_tot_no_ice = 3.0
        self.Q_tot_ice = 0.5
        self.T_ref = -1.8
        self.T_bot = -1.8
        self.B_t = 5.0
        self.n_w = 2.5
        self.n_s = 2.8
        self.k_c = 2.0
        self.C_0 = 20.0
        self.k_t = np.array([[0, 10**(-4),   10**(-5), 10**(-7)],    #ICE     ice-ML, ML-PC, PC-DP, DP-AB
                             [0, 6*10**(-4), 10**(-5), 10**(-7)]])   #NO ICE  ice-ML, ML-PC, PC-DP, DP-AB 
                            
        self.k_t_horiz = np.array([[0, 10**(-4),   10**(-5), 10**(-7)],    #ICE     ice-ML, ML-PC, PC-DP, DP-AB
                             [0, 6*10**(-4), 10**(-5), 10**(-7)]])   #NO ICE  ice-ML, ML-PC, PC-DP, DP-AB 
                             
#        self.k_t_horiz =1 * np.array([[10**(4), 10**(4),   10**(5), 10**(7)],    #ICE     ice-ML, ML-PC, PC-DP, DP-AB
#                                   [6*10**(4), 6*10**(4), 10**(5), 10**(7)]])   #NO ICE  ice-ML, ML-PC, PC-DP, DP-AB 
        
        self.L_f = 297000 #[J/kg] #heat of fusion for seawater. Taken as constant S=4 [0/00] T = -1.8 [C]

        #grid dependent UPDATE
        self.k_s_tot = np.ones(shape=(dim.n_x, dim.n_y, self.layers)) * 1 #Update only
        self.k_t_tot = np.ones(shape=(dim.n_x, dim.n_y, self.layers)) * 1 #update only
        
        self.k_s_horiz_tot = np.ones(shape=(dim.n_x, dim.n_y, self.layers)) * 1 #Update only
        self.k_t_horiz_tot = np.ones(shape=(dim.n_x, dim.n_y, self.layers)) * 1 #update only
        
        self.W_brine = np.zeros(shape=(dim.n_x, dim.n_y))
        self.W_melt = np.zeros(shape=(dim.n_x, dim.n_y)) 
        self.W_int = np.zeros(shape=(dim.n_x, dim.n_y, self.layers))  
        self.W_horiz = np.zeros(shape=(dim.n_x, dim.n_y, self.layers))
        
        self.Q_tot = np.ones(shape=(dim.n_x, dim.n_y, self.layers))
        self.Q_trans = np.zeros(shape=(dim.n_x, dim.n_y, layers)) #[W/m^2]        
        self.Q_int = np.zeros(shape=(dim.n_x, dim.n_y, layers))
        self.Q_turb = np.zeros(shape=(dim.n_x, dim.n_y))
        self.Q_horiz = np.zeros(shape=(dim.n_x, dim.n_y, layers))

   
        #FUNCTIONS      
        #ODE
        
        self.temp_func = temp_func
        self.salt_func = salt_func
        #salinity
        self.W_brine_func = W_brine_func
        self.W_melt_func = W_melt_func
        self.W_int_func = W_int_func      
        self.W_export_func = W_export_func
        self.W_tot_func = W_tot_func
        self.W_horiz_func = W_horiz_func
        #temperature
        self.Q_trans_func = Q_trans_func
        self.Q_int_func = Q_int_func
        self.Q_turb_func = Q_turb_func
        self.Q_tot_func = Q_tot_func
        self.Q_horiz_func = Q_horiz_func
       
        self.p_func = p_func
        self.c_w_func = c_w_func
    
        
        
        #update
        self.update_func = update_func

  
    def new_value(self,
                  start = 0,
                  **kw):
                      
        return np.ones(shape=(self.dim.n_x, self.dim.n_y, self.layers)) * start
#      
    def new_value_2(self,
                   name,
                   start = 0,
                   **kw):      
        self.name = np.ones(shape=(self.dim.n_x, self.dim.n_y, self.layers)) * start
        return self
    
class ice(object):
    '''Ice layer
    
    Attributes
    ----------
        C_r : float
            Ice export fraction
            
        p_ice : float
            Desnity of ice #[kg/m^3]
  
        V : numpy array, shape=(dim.n_x, dim.n_y)
            Ice volume
        
        h_ice : numpy array, shape=(dim.n_x, dim.n_y)
            Height of the ice-layer
        
        x_ice : numpy array, shape=(dim.n_x, dim.n_y)
            Fraction of the ocean covered in ice.             
            
        T_top : numpy array, shape=(dim.n_x, dim.n_y)
            Temperature on the top of the ice
            
        I_export : numpy array, shape=(dim.n_x, dim.n_y)
            Ice export rate.
        
        dV : numpy array, shape=(dim.n_x, dim.n_y)
            Ice grow rate (without I_export)
            
        export = numpy array, shape=(int(1/dim.dt_y), dim.n_x, dim.n_y))
            Ice grow rate for every step last year.
            
        dp_ice : numpy array, shape=(dim.n_x, dim.n_y)
            polynya ice
        
        Q_sw : numpy array, shape=(dim.n_x, dim.n_y)
            Absorbed long short wave radiation        
        
        Q_lw : numpy array, shape=(dim.n_x, dim.n_y)        
            Reflected long wave radiation            
        
        V_save : numpy array shape=(dim.steps+1, dim.n_x, dim.n_y))     
            save array
            
    Function pointers
    ----------------
        update_func : <type 'function'>
            Pointer to a function that updates all ice variables. 
        
        dp_ice_func : <type 'function'>
            Polynya ice.
            
        T_top_func : <type 'function'>
            Temperature on top of the ice.
        
        V_func : <type 'function'>
            Ice volume.
        
        dV_func : <type 'function'>
            Ice grow rate (without I_export).
            
        
    '''
    def __init__(self,
                 dim = dim(),
                 V_start = np.array([1.0]),
                 h_start = 4.0,
                 x_start = 0.5,
                 T_top_start = 0.0,
                 I_export_start = 0.0,
                 V_growth_start = 0.000,
                 C_r = 0.05,
                 update_func = update.ice_singh,
                 V_func = ice_volume.singh,
                 dp_ice_func = temperature.dp_ice,
                 T_top_func = temperature.T_top,
                 dV_func = ice_volume.dV,
                 **kw):
   
        #constants
        self.C_r = C_r   
        self.p_ice = 910.0 #[kg/m^3] density ice SHOULD MAYBE BE CHECKED AGAIN! 


        #grid dependent variables       
        self.V = np.ones(shape=(dim.n_x, dim.n_y)) * V_start
        self.h_ice = np.ones(shape=(dim.n_x, dim.n_y)) * h_start
        self.x_ice = np.ones(shape=(dim.n_x, dim.n_y)) * x_start
        self.T_top = np.ones(shape=(dim.n_x, dim.n_y)) * T_top_start
        self.I_export = np.ones(shape=(dim.n_x, dim.n_y)) * I_export_start #NOT CORRECT
        self.dV = np.ones(shape=(dim.n_x, dim.n_y)) * V_growth_start #growth/melt rate (without export) #SHOULD IMPORT BE PART??????????????????????????
        self.export = np.ones(shape=(int(1/dim.dt_y), dim.n_x, dim.n_y)) * V_growth_start #
        self.dp_ice = np.zeros(shape=(dim.n_x, dim.n_y))
        self.Q_sw = np.zeros(shape=(dim.n_x, dim.n_y))
        self.Q_lw = np.zeros(shape=(dim.n_x, dim.n_y))          
            
        #year dependent
        self.V_save = np.zeros(shape=(dim.save_steps+1, dim.n_x, dim.n_y)); self.V_save[0,:,:] = self.V         
        self.h_ice_save = np.zeros(shape=(dim.save_steps+1, dim.n_x, dim.n_y)); self.h_ice_save[0,:,:] = self.h_ice         
        #Functions
        self.update_func = update_func
        self.dp_ice_func = dp_ice_func
        self.T_top_func = T_top_func
        self.V_func = V_func
        self.dV_func = dV_func
    
if __name__ == "__main__":
    '''
    Only for testing
    '''
#    ode.singh()
    
#    dim_test = dim()
#    world_test = world(dim = dim_test)
  
#    world_test.plot_all()
#    world_test.animation(world_test.ocean.T_save)    
#    print world_test.ice.export
    
        
