# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 09:44:40 2016

@author: irg16
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

x_initial = 7.*R_m
y_initial = 3.*R_m

escape_v = ((2*G*M_m)/(((x_initial)**2) + ((y_initial)**2))**0.5)**0.5
esc_v_split = escape_v / (2**0.5)
print escape_v
print esc_v_split 

for k in np.arange(0, 2000, 500):
    initial=[x_initial, y_initial, -k, k]
    soln=spi.odeint(f,initial,t)  
    x=soln[:,0]
    y=soln[:,1]
    vx=soln[:,2]
    vy=soln[:,3]
    
    R=((x**2)+(y**2))**0.5 
    Vsq=((vx**2)+(vy**2))
    
    PE = -G*M_m*M_sat/(R)
    KE = 0.5*M_sat*Vsq

    pl.figure(1)
    pl.plot(t, PE)
    pl.xlabel("Time (s)")
    pl.ylabel("Potential Energy")
    
    pl.figure(2)
    pl.plot(t, KE)
    pl.xlabel("Time (s)")
    pl.ylabel("Kinetic Energy")
    
    pl.figure(3)
    pl.plot(t, KE + PE)
    pl.xlabel("Time (s)")
    pl.ylabel("Total Energy")

    pl.figure(4)    
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