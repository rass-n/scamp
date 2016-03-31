# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:13:23 2016

@author: ko
"""

import scamp as sc
import numpy as np
import matplotlib.pyplot as plt
import time

#dim = sc.dim()
#print 'steps', dim.steps
#ocean = sc.ocean(dim, ode_func=sc.ode.test)
#ocean.test = ocean.new_value(4)
#ocean.ode_func = sc.ode.singh
#world = sc.world(dim, ocean_start=ocean)
#
#
#
#print world.ocean.ode_func
#sc.forward_euler(world)
#print world.ocean.T

np.array	
#simple test
#dim = sc.dim(dt = 1/2.0)
dim = sc.dim(n_x = 1, n_y = 1, years = 100, dt=1/20.0)
#ocean = sc.ocean(dim, temp_func = sc.temperature.singh_gl, update_func=sc.update.ocean_singh_gl, ode_func = sc.ode.singh_gl)
ocean = sc.ocean(dim)
atmos = sc.atmos(dim)
ice = sc.ice(dim, C_r = 0.8)
#ice.test = 0

world = sc.world(dim, ocean_start = ocean, ice_start = ice, atmos_start = atmos)
print type(world.ocean.Q_turb_func)

#world = world.update_func(world, 1)

#print world.ocean.m / world.ocean.c_w



world = sc.ode_solver.forward_euler(world)



#world.animation(world.ocean.T_save)
#print world.ice.test

#world.plot_con()

#sum_year = np.zeros(world.dim.years)
##print 'HER', world.dim.years
#for i in xrange(world.dim.years):
#    for k in xrange(20):
##        print 20.0*i + k
#        sum_year[i] += sum(world.ocean.S_save[20*i+k,0,0,:] * world.ocean.p_save[20*i+4,0,0,:] * world.ocean.h)
#
#
#plt.plot(sum_year)
#plt.show()


    
#world.plot_all()
#print ' '
#print 'NY'
#print ' '

#dim1 = sc.dim()
#dim1 = sc.dim(dt = 1/2.0)
#world1 = sc.world(dim1)
#world1 = sc.ode_solver.forward_euler(world1)




#np.set_printoptions(precision=40)
#print world.ocean.T
#print world.ocean.S
#print world.ice.V
 


#T_test = np.array([[[-1.4599290386061647417648146074498072266579,
#    1.7022975696284661228929735443671233952045,
#    1.8659021856765445335213371436111629009247,
#    1.0000008608371402374359604436904191970825]]])
#    
#S_test = np.array([[[ 33.3480076248373791258927667513489723205566,
#    33.5432154633661241405206965282559394836426,
#    33.4856470750004007186362287029623985290527,
#    34.9999998947139019378482771571725606918335]]])
#    
#V_test = np.array([[ 2.4731354557928586501702739042229950428009]])
#
#print 'TEST'
#print T_test - world.ocean.T
#print S_test - world.ocean.S
#print V_test - world.ice.V



##print world.ice.dV
##print "START"
##world.print_all()
##
#time1 = time.time()

##


#dim1 = sc.dim(dt = 1/2.0)
#dim1 = sc.dim()
##ocean1 = sc.ocean(dim1, salt_func = sc.salinity.singh_gl)
#ocean1 = sc.ocean(dim1)
#atmos1 = sc.atmos(dim1)
#ice1 = sc.ice(dim1)
#world1 = sc.world(dim1, ocean_start = ocean1, ice_start = ice1, atmos_start = atmos1)
#print  ' '
#print 'NY AMOK'
#world1 = sc.ode_solver.forward_euler(world1)
#
#
#print 'ER DE ENS'
#print world1.ocean.T - world.ocean.T
#print world1.ocean.S - world.ocean.S
#print world1.ice.V - world.ice.V



#time2 = time.time()
##print "\nSLUT"
##world.print_all()
#world.plot_all()

#world.print_all()

#T_test = np.array([[[-1.4599309692026485318905315580195747315883636474609375,
#    1.702298760372614072622354797204025089740753173828125,
#    1.8659021835176721193505500195897184312343597412109375,
#    1.0000008608371431240158244690974242985248565673828125]]])
#    
#S_test = np.array([[[ 33.313715678019207189208827912807464599609375,
#    33.538125303658461007216828875243663787841796875,
#    33.4856518150587589843780733644962310791015625,
#    34.99999989471398720297656836919486522674560546875]]])
#V_test = np.array([[ 2.47313708799727915987887172377668321132659912109375]])







##
#import sys
#sys.path.append('/home/ko/Master/backward/Runge-kutta')    
#import forward_euler
#import numpy
###
#dim2 = sc.dim()
#world2 = sc.world(dim2)
#
###
#time3 = time.time()
#[T_t, S_t, V_t, p_t] = forward_euler.forward_euler_test(world2.ocean.T, world2.ocean.S, world2.ice.V, world2.ocean.h, world2.dim.dt_y, world2.dim.years)
#time4 = time.time()
#
#
#
#print 'Tid min', time2-time1
#
#print 'tid gam', time4-time3
#time = numpy.arange(len(S_t[:,0,0,0])) * world2.dim.dt_y
#
#plt.plot(time,S_t[:,0,0,0], 'r', label='ML'), plt.plot(time,S_t[:,0,0,1],'b', label='PC') , plt.plot(time,S_t[:,0,0,2],'c', label='DP')
#plt.plot(time,S_t[:,0,0,3],'g', label='AB')
#plt.xlabel('time (steps)')
#plt.ylabel('Salinity')
#plt.title('Salinity Layers')
#plt.grid(True)
#legend = plt.legend(loc='upper left', shadow=True, fontsize='x-large')
#plt.show()
#
#plt.plot(time,T_t[:,0,0,0], 'r', label='ML'), plt.plot(time,T_t[:,0,0,1],'b', label='PC') , plt.plot(time,T_t[:,0,0,2],'c', label='DP')
#plt.plot(time,T_t[:,0,0,3],'g', label='AB')
#plt.xlabel('time (steps)')
#plt.ylabel('Temperature')
#plt.title('Temperature')
#plt.grid(True)
#legend = plt.legend(loc='upper left', shadow=True, fontsize='x-large')
#plt.show()
#
#    
#plt.plot(time,V_t[:,0,0], 'r', label='ML'), plt.plot(time,V_t[:,0,0],'b', label='PC') , plt.plot(time,V_t[:,0,0],'c', label='DP')
#plt.plot(time,V_t[:,0,0],'g', label='AB')
#plt.xlabel('time (steps)')
#plt.ylabel('ice_vol')
#plt.title('ice_volume')
#plt.grid(True)
#legend = plt.legend(loc='upper left', shadow=True, fontsize='x-large')
#plt.show()
#
#
#plt.plot(time,p_t[:,0,0,0], 'r', label='ML'), plt.plot(time,p_t[:,0,0,1],'b', label='PC') , plt.plot(time,p_t[:,0,0,2],'c', label='DP')
#plt.plot(time,p_t[:,0,0,3],'g', label='AB')
#plt.xlabel('time (steps)')
#plt.ylabel('Density')
#plt.title('Density Layers')
#plt.grid(True)
#legend = plt.legend(loc='upper left', shadow=True, fontsize='x-large')
#plt.show()


