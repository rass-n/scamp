# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 18:15:53 2016

@author: ko
"""

import matplotlib.animation as animation
import matplotlib.pyplot as plt


def ani(X, steps):
    print 'steps', steps
    fig = plt.figure()
#
#
##
    im1 = plt.imshow(X[0,:,:,1], vmin=min(X[:,0,0,1]), vmax=max(X[:,0,0,1]))
    plt.colorbar()
    def updatefig1(i):
        print 'i', i
        im1.set_data(X[i+1,:,:,1])
        return im1,

    ani = animation.FuncAnimation(fig, updatefig1, frames=steps, blit=True,repeat=False)
    plt.show()
    
    return
# 
#im2 = plt.imshow(V[0,:,:], vmin=min(V[:,0,0]), vmax=max(V[:,0,0]))
#plt.colorbar()
#def updatefig2(i):
#    im2.set_data(V[i+1,:,:])
#    return im2,
#
#ani = animation.FuncAnimation(fig, updatefig2, frames=steps, blit=True,repeat=False)
#plt.show()

 
#im3 = plt.imshow(S[0,5,:,:], vmin=30, vmax=40)
#plt.colorbar()
#def updatefig3(i):
##    print i
#    im3.set_data(S[i+1,5,:,:])
#    return im3,
##
#ani = animation.FuncAnimation(fig, updatefig3, frames=steps, blit=True,repeat=False)
#plt.show()
###
    