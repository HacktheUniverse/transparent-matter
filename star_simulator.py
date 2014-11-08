# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 21:55:12 2014

@author: lucio tolentino

Script for calculating the movement of stars under different amount of dark
matter.


"""

import sys
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
import json

#Will be called with input as JSON
parameters = sys.argv[0]
parameters = {"dark_matter":True, "amount_dark_matter": 6, "distribution": "Scenario A"}

width = 50
height = 50
number_of_stars = 7

#place the stars randomly
stars_x =  [random.uniform(width) for i in range(number_of_stars)] # [1, 0.5, 1.5] #
stars_y = [random.uniform(height) for i in range(number_of_stars)] # [0, 0.5, 0.5] # 

#move them a little bit
timesteps = 128
theta = np.pi/64
positions = {0:{i:(stars_x[i], stars_y[i]) for i in range(number_of_stars)}} # save initial positions
for t in range(1,timesteps):
    positions[t] = {}
    for s in range(number_of_stars):
        #rotate star position
        x = stars_x[s]
        y = stars_y[s]
        stars_x[s] = x*np.cos(theta) - y*np.sin(theta)
        stars_y[s] = y*np.cos(theta) + x*np.sin(theta)
        
        #save position
        positions[t][s] = (stars_x[s], stars_y[s])

#visualize
colors = ['r', 'b', 'g', 'k', 'y', 'c', 'm',]
for s in range(number_of_stars):
    x = [positions[t][s][0] for t in range(timesteps)]
    y = [positions[t][s][1] for t in range(timesteps)]
    #print "x", x
    #print "y", y
    plt.scatter(x,y,c = colors[s], linewidth=0)
    plt.xlim((-width*1.5,width*1.5))
    plt.ylim((-height*1.5,height*1.5))
#print positions
    
json.dumps(positions)
    
    