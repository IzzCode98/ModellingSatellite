# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

G=6.67*(10**-11) #m^3 kg^-1 s^-2 #gravitational constant
M_m=6.4*(10**23) #kg #mass of mars
R_m=3.4*(10**6) #m #radius of mars
M_sat=260 #kg #mass of satellite

import scipy as sp
import numpy as np
import pylab as pl
import scipy.integrate as spi

def f(value,t):
    x_x=value[0]
    x_y=value[1]
    v_x=value[2]
    v_y=value[3]    
    Rsq=(x_x**2)+(x_y**2)    
    a_x=-G*M_m*x_x/(Rsq)**1.5
    a_y=-G*M_m*x_y/(Rsq)**1.5
    if Rsq > R_m**2:
        return [v_x, v_y, a_x, a_y] 
    else:
        return [0, 0, 0, 0]
    

t=sp.linspace(0.,1000000.,1000) 

#create k for loop
# set initial conditions
for k in np.arange(0, 10000, 100):
    initial=[7.*R_m, 3.*R_m, -k, k]
    soln=spi.odeint(f,initial,t) 
        
    x=soln[:,0]
    y=soln[:,1]
    vx=soln[:,2]
    vy=soln[:,3]
    
    #pl.figure(1)
    #pl.plot(t, x)
    #pl.xlabel("time (s)")
    #pl.ylabel("x displacement (m)")
    #
    #pl.figure(2)
    #pl.plot(t, y)
    #pl.xlabel("time (s)")
    #pl.ylabel("y displacement (m)")
    #
    #pl.figure(3)
    #pl.plot(t, vx)
    #pl.xlabel("time (s)")
    #pl.ylabel("x velocity (m/s)")
    #
    #pl.figure(4)
    #pl.plot(t, vy)
    #pl.xlabel("time (s)")
    #pl.ylabel("y velocity (m/s)")
    
    #pl.figure(5)
    pl.plot(x, y)
    pl.xlabel("x displacement (m)")
    pl.ylabel("y displacement (m)")
    pl.axis([-0.5e8,0.5e8,-0.5e8,0.5e8])
    def c(q, r):
        y = ((r**2)-(q**2))**0.5    
        return y  
    def cm(q, r):
        y = -((r**2)-(q**2))**0.5    
        return y    
    q = np.arange(-R_m, R_m, 100)
    pl.plot(q, c(q, R_m), 'r-')
    pl.plot(q, cm(q, R_m), 'r-')
    
pl.show()



#escape_v = ((G*M_m)/((7.*R_m)**2 + (3.*R_m)**2)**0.5)**0.5
#ang = np.arctan(3./7.)
#evx = escape_v * np.cos(ang)
#evy = escape_v * np.sin(ang)
#print escape_v, evx, evy       
#
#initial=[7.*R_m, 3.*R_m, -evx, evy]:
#soln=spi.odeint(f,initial,t) 
#    
#x=soln[:,0]
#y=soln[:,1]
#vx=soln[:,2]
#vy=soln[:,3]
#
#pl.plot(x, y)
#pl.xlabel("x displacement (m)")
#pl.ylabel("y displacement (m)")
#pl.axis([-0.5e8,0.5e8,-0.5e8,0.5e8])
#def c(q, r):
#    y = ((r**2)-(q**2))**0.5    
#    return y  
#def cm(q, r):
#    y = -((r**2)-(q**2))**0.5    
#    return y    
#q = np.arange(-R_m, R_m, 100)
#pl.plot(q, c(q, R_m), 'r-')
#pl.plot(q, cm(q, R_m), 'r-')
#
#pl.show()

#when it isn’t captured, 
    #what is the angular deviation of the satellite’s trajectory 
    #and the distance of closest approach 
    #(you should plot the deviation as a function of initial velocity)