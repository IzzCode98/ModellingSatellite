# -*- coding: utf-8 -*-
# warm up task for Project C
"""
Spyder Editor

This is a temporary script file.
"""
#WARM UP

#SciPy has several functions for solving sets of first order ODE given a set of initial conditions. 
#The simplest of these (and the one we recommend you use) is odeint, which is part of the scipy.integrate package.

#Many physical problems are described by higher (at least second) order ODE, 
#however these can always be written down as a set of first order ODE. 
#For example, through Newton’s second law, 
#the position of an object is related to the force that it experiences by a second order ODE, 
#however it can be written as a set of two first order ODE. I.e. for a force F(x) on a body of mass m we have
    #d2x/dt2=F(x)m
#but this can be rewritten as the set of first order ODE
    #dvdt=F(x)m
    #dxdt=v

#In order to use the ODE integrator you will have to define a function whose calling arguments are variables to be differentiated 
#and the variable with respect to which they will be differentiated. 
#You will also need to specify the initial conditions.

#EXAMPLE CODE
#Let us consider the simple example of a particle moving in one dimension under gravity at the earth’s surface,
# with an initial velocity of +50ms^−1 
#(which we define to be upwards so that gravity is in the negative direction). 
#We might simulate this with the following code:

import scipy as sp
import numpy as np
import pylab as pl
import scipy.integrate as spi
#defining imported modules with abbreviated names to make life easier

g=9.8 #define local gravity

#called with the variables to be differentiated
def f(input_data,t): #this defines a function f with parameters input_data and t to have the following properties
  loc=input_data[0] #input_data is a list and we can define the values in the list by assigning them, e.g. the first value is loc
  speed=input_data[1] #the second value in the input_data list is speed
  acc=-g #don't forget it is in the negative direction #acc is simply a constant in this scenario
  return [speed,acc]
#returns differentiated values
  
  #see notes about defining functions in google Doc


t=sp.linspace(0.,10.,100) # solving every tenth of a second - more than is needed 
            #the last number (100) represents the number of data points within the given range
            #this defines t for the soln equation

# set initial conditions going up at 50 ms^-1 starting at ground level

initial_values=[0.,50] 
            #this defines initial_values for the soln equation

soln=spi.odeint(f,initial_values,t) #this odeint thingy is arranged such that it looks at function f and begins the integration process with the initial_values at steps of t - all of which have already been defined


print soln #soln now has 2 columns, the first being the values of loc for every t step and the second being the values of speed for every t step - these were generated by the above equation using the f loop


x=soln[:,0] #this assigns the first column of soln (including all rows via :) to a value named x
v=soln[:,1] #this assigns the second column of soln (including all rows via :) to a value named v

pl.figure(1)
pl.plot(t,x)
pl.xlabel("time (s)")
pl.ylabel("height (m)")

pl.figure(2)
pl.plot(t,v)
pl.xlabel("time (s)")
pl.ylabel("velocity (m/s)")


pl.show() #now this information has been plotted