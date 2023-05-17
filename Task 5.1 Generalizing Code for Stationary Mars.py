# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 11:38:50 2016

@author: irg16
"""

#Generalizing all code. Task 5.


G=6.67*(10**-11) #m^3 kg^-1 s^-2 #gravitational constant
M_m=6.4*(10**23) #kg #mass of mars
R_m=3.4*(10**6) #m #radius of mars
M_sat=260 #kg #mass of satellite

import scipy as sp
import numpy as np
import pylab as pl
import scipy.integrate as spi

pl.close(pl.figure(1))
pl.close(pl.figure(2))
pl.close(pl.figure(3))
pl.close(pl.figure(4))
pl.close(pl.figure(5))
pl.close(pl.figure(6))
pl.close(pl.figure(7))


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


t=sp.linspace(0.,1000000.,10000) 

example_velocity = 1000

x_initial = 7.*R_m
y_initial = 3.*R_m

escape_v = ((2*G*M_m)/(((x_initial)**2) + ((y_initial)**2))**0.5)**0.5
esc_v_split = escape_v / (2**0.5) #because x and y velocities are of the same magnitude

dist = []
collide = []
v_list = []
min_list = []

del min_list [:]
del v_list [:]
del dist [:]
del collide [:]


for k in np.arange(0, 10000, 100):
    
    initial=[x_initial, y_initial, -k, k]
   
    soln=spi.odeint(f,initial,t) 
    
     
    if k == example_velocity:
      
        x=soln[:,0]
        y=soln[:,1]
        vx=soln[:,2]
        vy=soln[:,3]
        
        pl.figure(1)
        pl.subplot(3,2,1)
        pl.plot(t, x, 'b')
        pl.xlabel("time (s)")
        pl.ylabel("x displacement (m)")
        
        pl.subplot(3,2,2)
        pl.plot(t, y, 'b')
        pl.xlabel("time (s)")
        pl.ylabel("y displacement (m)")
        
        pl.subplot(3,2,3)
        pl.plot(t, vx, 'm')
        pl.xlabel("time (s)")
        pl.ylabel("x velocity (m/s)")
        
        pl.subplot(3,2,4)
        pl.plot(t, vy, 'm')
        pl.xlabel("time (s)")
        pl.ylabel("y velocity (m/s)")
        
        pl.subplot(3,2,5)
        pl.plot(x, y, 'g')
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
        
    else:
        x=soln[:,0]
        y=soln[:,1]
        vx=soln[:,2]
        vy=soln[:,3]
        
        pl.figure(2)
        pl.subplot (2,2,1)
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

        if k <= esc_v_split:
                        
            x=soln[:,0]
            y=soln[:,1]
            vx=soln[:,2]
            vy=soln[:,3]
            
            pl.figure(2)
            pl.subplot (2,2,3)            
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
                        
            del collide [:]
        
            collide = [(((a**2)+(b**2))**0.5) for a, b in zip(x, y)]
            
            vcollide = k * (2**0.5)  
                        
            if min(collide) <= R_m:
                print 'Satellite has crashed: it had initial velocity of', vcollide 
        
        if k >= esc_v_split:
                          
            x=soln[:,0]
            y=soln[:,1]
            vx=soln[:,2]
            vy=soln[:,3]
            

            pl.figure(2)
            pl.subplot (2,2,4)            
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
            
            
            del dist[:]

            dist = [(((a**2)+(b**2))**0.5) for a, b in zip(x, y)]    
            

                
            min_list.append (min(dist))

            vinit= k * (2**0.5)  
            v_list.append(vinit)


        
            vxfinal=soln[99,2]
            vyfinal=soln[99,3]
            angfinal=np.arctan(vyfinal/vxfinal)
            anginitial=np.arctan(1)
            angledev=anginitial+angfinal
            
            pl.figure(3)
            pl.subplot(2,2,1)
            pl.plot (vinit, angledev , 'r+', label='Angle of Deviation')
            pl.ylabel('Angle of Deviation (rad)')
            pl.xlabel('Initial Velocity (m/s)')
            pl.axvline (x=escape_v, color='b', linestyle='--', label='Escape Velocity')
        
            pl.figure(3)
            pl.subplot(2,2,3)    
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
            

pl.figure(3)   
pl.subplot(2,2,2)
pl.plot (v_list, min_list , 'r-', label='Distance of Closest Approach')
pl.axvline (x=escape_v, color='b', linestyle='--', label='Escape Velocity')
pl.ylabel('Distance of Closest Approach (m)')
pl.xlabel('Initial Velocity (m/s)')
pl.legend(loc='lower right')

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

    pl.figure(4)
    pl.subplot(2,2,1)
    pl.plot(t, PE)
    pl.xlabel("Time (s)")
    pl.ylabel("Potential Energy")
    
    pl.figure(4)
    pl.subplot(2,2,2)
    pl.plot(t, KE)
    pl.xlabel("Time (s)")
    pl.ylabel("Kinetic Energy")
    
    pl.figure(4)
    pl.subplot(2,2,3)
    pl.plot(t, KE + PE)
    pl.xlabel("Time (s)")
    pl.ylabel("Total Energy")

    pl.figure(4)
    pl.subplot(2,2,4)
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

#Moving Mars

def f2(value,t):
    x_x=value[0]
    x_y=value[1]
    v_x=value[2]
    v_y=value[3]    
    x_mars= 24100 * t
    y_mars= 0
    Rsq=((x_x-x_mars)**2)+((x_y-y_mars)**2)    
    a_x=-G*M_m*(x_x-x_mars)/(Rsq)**1.5
    a_y=-G*M_m*(x_y-y_mars)/(Rsq)**1.5
    if Rsq > R_m**2:
        return [v_x, v_y, a_x, a_y] 
    else:
        return [0, 0, 0, 0]

x_mars= 24100 * t

for k in np.arange(0, 2000, 100):
    
    initial=[x_initial, y_initial, 24100-k, k]
   
    soln=spi.odeint(f2,initial,t) 
    
     
    if k == 100:
      
        x=soln[:,0]
        y=soln[:,1]
        vx=soln[:,2]
        vy=soln[:,3]
        
        pl.figure(5)
        pl.subplot(3,2,1)
        pl.plot(t, x, 'b')
        pl.xlabel("time (s)")
        pl.ylabel("x displacement (m)")
        
        pl.subplot(3,2,2)
        pl.plot(t, y, 'b')
        pl.xlabel("time (s)")
        pl.ylabel("y displacement (m)")
        
        pl.subplot(3,2,3)
        pl.plot(t, vx, 'm')
        pl.xlabel("time (s)")
        pl.ylabel("x velocity (m/s)")
        
        pl.subplot(3,2,4)
        pl.plot(t, vy, 'm')
        pl.xlabel("time (s)")
        pl.ylabel("y velocity (m/s)")
        
        pl.subplot(3,2,5)
        pl.plot(x, y, 'g')
        pl.xlabel("x displacement (m)")
        pl.ylabel("y displacement (m)")
        pl.axis([0e8,2e10,-0.5e8,0.5e8])
        
    if k == 1500:
      
        x=soln[:,0]
        y=soln[:,1]
        vx=soln[:,2]
        vy=soln[:,3]
        
        pl.figure(6)
        pl.subplot(3,2,1)
        pl.plot(t, x, 'b')
        pl.xlabel("time (s)")
        pl.ylabel("x displacement (m)")
        
        pl.subplot(3,2,2)
        pl.plot(t, y, 'b')
        pl.xlabel("time (s)")
        pl.ylabel("y displacement (m)")
        
        pl.subplot(3,2,3)
        pl.plot(t, vx, 'm')
        pl.xlabel("time (s)")
        pl.ylabel("x velocity (m/s)")
        
        pl.subplot(3,2,4)
        pl.plot(t, vy, 'm')
        pl.xlabel("time (s)")
        pl.ylabel("y velocity (m/s)")
        
        pl.subplot(3,2,5)
        pl.plot(x, y, 'g')
        pl.xlabel("x displacement (m)")
        pl.ylabel("y displacement (m)")
        pl.axis([0e8,2e10,-0.5e8,0.5e8])
        
    if k == 1000:
      
        x=soln[:,0]
        y=soln[:,1]
        vx=soln[:,2]
        vy=soln[:,3]
        
        pl.figure(7)
        pl.subplot(3,2,1)
        pl.plot(t, x, 'b')
        pl.xlabel("time (s)")
        pl.ylabel("x displacement (m)")
        
        pl.subplot(3,2,2)
        pl.plot(t, y, 'b')
        pl.xlabel("time (s)")
        pl.ylabel("y displacement (m)")
        
        pl.subplot(3,2,3)
        pl.plot(t, vx, 'm')
        pl.xlabel("time (s)")
        pl.ylabel("x velocity (m/s)")
        
        pl.subplot(3,2,4)
        pl.plot(t, vy, 'm')
        pl.xlabel("time (s)")
        pl.ylabel("y velocity (m/s)")
        
        pl.subplot(3,2,5)
        pl.plot(x, y, 'g')
        pl.xlabel("x displacement (m)")
        pl.ylabel("y displacement (m)")
        pl.axis([0e8,2e10,-0.5e8,0.5e8])
        
        
for k in np.arange(100, 2000, 400):
    initial=[x_initial, y_initial, 24100-k, k]
    soln=spi.odeint(f2,initial,t)  
    x=soln[:,0]
    y=soln[:,1]
    vx=soln[:,2]
    vy=soln[:,3]
    
    R=(((x-x_mars)**2)+(y**2))**0.5 
    Vsq=((vx**2)+(vy**2))
    
    PE = -G*M_m*M_sat/(R)
    KE = 0.5*M_sat*Vsq

    pl.figure(8)
    pl.subplot(2,2,1)
    pl.plot(t, PE)
    pl.xlabel("Time (s)")
    pl.ylabel("Potential Energy")
    
    pl.subplot(2,2,2)
    pl.plot(t, KE)
    pl.xlabel("Time (s)")
    pl.ylabel("Kinetic Energy")
    
    pl.subplot(2,2,3)
    pl.plot(t, KE + PE)
    pl.xlabel("Time (s)")
    pl.ylabel("Total Energy")

    pl.subplot(2,2,4)
    pl.plot(x, y)
    pl.xlabel("x displacement (m)")
    pl.ylabel("y displacement (m)")
    pl.axis([0e8,2e10,-0.5e8,0.5e8])
    
pl.show()