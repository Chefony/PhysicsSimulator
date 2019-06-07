# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:48:47 2019

@author: Stefan Ivanov
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import timeit

dx = 0.1
dy = 0.1
start = timeit.default_timer()
def make_hot_corner_plate(dx=0.01,dy=0.01):
# =============================================================================
#     Creates a 2D araay representin a plate with corner
#    temperature of 30C and center temperature of 0C
# =============================================================================
    x = np.arange(-1,1,dx)
    y = np.arange(-1,1,dy) 
    x, y = np.meshgrid(x,y)
    return (x**2+y**2)*15

def get_derivative_mirror_boundary(T,x,y,a,h):
# =============================================================================
#     Returns the second derivative with respect to space
#     T-2D array containing temperatures
#     x,y-coordinates
#     a-alpha(heat diffusivity)
#     h-step size
# =============================================================================
    dTdx = 0.0
    if(x==0):
        if(y==0):
            dTdx = (T[x+1,y]+T[x,y+1]-2*T[x,y])/(2*h**2)
        elif(y==T.shape[1]-1):
            dTdx = (T[x+1,y]+T[x,y-1]-2*T[x,y])/(2*h**2)
        else:
            dTdx = (2*T[x+1,y]+T[x,y-1]+T[x,y+1]-4*T[x,y])/h**2
    elif(x==T.shape[0]-1):
        if(y==0):
            dTdx = (T[x-1,y]+T[x,y+1]-2*T[x,y])/(2*h**2)
        elif(y==T.shape[1]-1):
            dTdx = (T[x-1,y]+T[x,y-1]-2*T[x,y])/(2*h**2)
        else:
            dTdx = (2*T[x-1,y]+T[x,y-1]+T[x,y+1]-4*T[x,y])/h**2
    elif(y==0):
        if(x==0 or x==T.shape[0]-1):
            pass
        else:
            dTdx = (T[x-1,y]+T[x+1,y]+2*T[x,y+1]-4*T[x,y])/h**2
    elif(y==T.shape[1]-1):
        if(x==0 or x==T.shape[0]-1):
            pass
        else:
            dTdx = (T[x-1,y]+T[x+1,y]+2*T[x,y-1]-4*T[x,y])/h**2
    else:
        dTdx = (T[x-1,y]+T[x+1,y]+T[x,y-1]+T[x,y+1]-4*T[x,y])/h**2
    return a*dTdx

def time_step(T,a,dt,h):
# =============================================================================
#     Returns a 2-D array with temperatures after a period of time(dt)
#     calculated from T using h as step size ana a as heat diffusivity
# =============================================================================
    result = np.empty((T.shape[0],T.shape[1]))
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            result[i,j] = T[i,j] + get_derivative_mirror_boundary(T,i,j,a,h)*dt
    return result

def time_interval(T,alpha,dt,epsilon,t_final):
# =============================================================================
#    Returns a set of 2-D arrays representing the change of temperatures(T)
# =============================================================================
    result = np.empty(int(t_final/dt),dtype=np.ndarray)
    result[0] = T
    for i in range(1,result.shape[0]):
        result[i] = time_step(result[i-1],alpha,dt,epsilon)
    return result
                         
dt = 0.005
time_steps = time_interval(make_hot_corner_plate(dx,dy),0.143,dt,dx,12)
fig = plt.figure()
ax = plt.gca(projection="3d")
X = np.arange(-1,1,dx)
Y = np.arange(-1,1,dy)
X, Y = np.meshgrid(X,Y)
Z = time_steps[0]
surf = ax.plot_surface(X, Y, time_steps[2399], cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

end = timeit.default_timer()

print(end-start)