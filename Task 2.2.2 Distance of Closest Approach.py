# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 16:01:11 2016

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
        
        
t=sp.linspace(0.,100000.,10000) 

x_initial = 7.*R_m
y_initial = 3.*R_m

escape_v = ((2*G*M_m)/(((x_initial)**2) + ((y_initial)**2))**0.5)**0.5
esc_v_split = escape_v / (2**0.5)
print escape_v
print esc_v_split

dist = []
v_list = []
min_list = []

del min_list [:]
del v_list [:]
del dist [:]


for k in np.arange(esc_v_split, 20000, 200):
    initial=[x_initial, y_initial, -k, k]
    soln=spi.odeint(f,initial,t)  
    #print soln
    x=soln[:,0]
    y=soln[:,1]
    vx=soln[:,2]
    vy=soln[:,3]
    
    del dist[:]

    dist = [(((a**2)+(b**2))**0.5) for a, b in zip(x, y)]    
    
    min_list.append (min(dist))
        #print dist
#    
#    vxinit=soln[0,2]
#    vyinit=soln[0,3]

    vinit= k * (2**0.5)  
    v_list.append(vinit)
    #print [vxinit, vyinit, vinit] 
    
#    vxfinal=soln[999,2]
#    vyfinal=soln[999,3]
#    angfinal=np.arctan(vyfinal/vxfinal)
#    anginitial=np.arctan(1)
#    angledev=anginitial+angfinal

    #print [angfinal, anginitial, angledev]
    #print [angledev , vinit]    
    
#    pl.figure(1)    
#    pl.plot (vinit, dist , 'r+')
#    pl.ylabel('Distance of Closest Approach (m)')
#    pl.xlabel('Initial Velocity (m/s)')
#    
    pl.figure(2)    
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

#print [v_list, min_list]

#print len (v_list)
#print len (min_list)
pl.figure(3)    
pl.plot (v_list, min_list , 'b-', label='Distance of Closest Approach')
pl.axvline (x=escape_v, color='r', linestyle='--', label='Escape Velocity')
pl.ylabel('Distance of Closest Approach (m)')
pl.xlabel('Initial Velocity (m/s)')
pl.legend(loc='lower right')



##to find closest apporach distance, 
