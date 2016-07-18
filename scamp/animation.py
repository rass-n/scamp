# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 18:15:53 2016

@author: Rasmus Nordfang
"""

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib


def ani(X, steps, vmin=30, vmax = 35, **kw):
    '''
    Simple animation of a full ocean-layer
    
    Input: X - numpy array with values. E.g. world.ocean.T 
           steps - int, number of frames     
    
    '''

    fig = plt.figure()
#
#
    my_cmap = matplotlib.cm.get_cmap('seismic')
    
#    im1 = plt.imshow(X[0,:,:,1], vmin=min(X[:,0,0,1]), vmax=max(X[:,0,0,1]), cmap = my_cmap)
    im1 = plt.imshow(X[0,:,:], vmin=vmin, vmax=vmax, cmap = my_cmap)
    plt.colorbar()
#    plt.colormaps('rainbow')
    def updatefig1(i):
#        print 'i', i
        im1.set_data(X[i+1,:,:])       
        return im1,
        
#    print 'Steps', steps
    ani = animation.FuncAnimation(fig, updatefig1, frames=steps, blit=True,repeat=False)
    plt.show()
    
    return
