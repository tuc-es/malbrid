#!/usr/bin/env python
# coding: utf-8

# # An object in an orbit around a planet -- damped version

# In[1]:


# Load Packages
import scipy
import matplotlib.pyplot
import numpy
import math
import malbrid
import random
random.seed(1234)


# In[2]:


def get_dynamics_and_zero_crossing_functions_orbit(state_name):
    AOrbit = numpy.array([[0,0,1,0,0],[0,0,0,1,0],[-1,0,-0.01,0,0],[0,-1,0,-0.01,0],[0,0,0,0,0]])
    
    xPos = simulator.get_var("xPos")
    yPos = simulator.get_var("yPos")
    xSpeed = simulator.get_var("xSpeed")
    ySpeed = simulator.get_var("ySpeed")
    const = simulator.get_var("const")                     
                         
    def nobump(x):
        return "OnlyOne",x,False
    
    if state_name=="OnlyOne":
        return AOrbit, []
    else:
        raise Exception("Internal Test error:"+str(state_name))

'''A test case for Randomized Dynamics -- Product state case'''
simulator = malbrid.LinearSystemSimulator(["xPos", "yPos", "xSpeed", "ySpeed", "const"])


# In[3]:


for i in range(0,1):
    simulator.simulate(get_dynamics_and_zero_crossing_functions_orbit, "OnlyOne",numpy.array([0,1,-1,0,1]),
                           max_time=1000,max_timestep=0.1)
    matplotlib.pyplot.plot(numpy.array(simulator.continuous_states)[:,0],numpy.array(simulator.continuous_states)[:,1])
    
            
# Finalize Plot
matplotlib.pyplot.xlabel('Time')
matplotlib.pyplot.ylabel('Temperature')
matplotlib.pyplot.show()


# In[4]:




