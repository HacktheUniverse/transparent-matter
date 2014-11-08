# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 21:55:12 2014

@author: lucio tolentino

Script for calculating the movement of stars under different amount of dark
matter.


"""

import sys
import os
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
import json
import scipy.constants as cons

#define movement functions -- given an x,y return x',y'
def rotate(x, y, theta):
    """
    Rotate a cartesian coordinate given by *x*, *y* by an angle *theta*.
    Returns the new cartesian coordinate *x'*, *y'*.
    """
    return x*np.cos(theta) - y*np.sin(theta), y*np.cos(theta) + x*np.sin(theta)

def angular_velocity(R, M):
    ly2m=1.0/1.05702341e-16
    return np.sqrt(cons.G*M/R/ly2m)

#Will be called with input as JSON
if __name__ == "__main__":
    #set parameters
    #parameters = sys.argv[0]
    parameters = {"dark_matter":False, "amount_dark_matter": 6, "distribution": "Scenario A"}
    width = 50
    number_of_stars = 20
    mass = 1e30
    timesteps = 100
    lambda_ = 0.2  # exponential distribution of stars
    
    #exponentially distribute stars from center of galaxy
    stars_r = [width*random.exponential(lambda_) for i in range(number_of_stars)]
    stars_theta = [random.uniform(2.0)*np.pi for i in range(number_of_stars)]
    stars_x = [stars_r[i]*np.cos(stars_theta[i]) for i in range(number_of_stars)]
    stars_y = [stars_r[i]*np.sin(stars_theta[i]) for i in range(number_of_stars)]
    
    #BUILD JSON
    positions = {0:{i:(stars_x[i], stars_y[i]) for i in range(number_of_stars)}} # save initial positions
    for t in range(1,timesteps):
        positions[t] = {}
        for s in range(number_of_stars):
            #rotate star's position
            if not parameters["dark_matter"]:
                new_angle = angular_velocity(stars_r[s], mass) * stars_r[s] * 2
                stars_x[s], stars_y[s] = rotate(stars_x[s], stars_y[s], new_angle)
            else:
                pass
            #save position
            positions[t][s] = (stars_x[s], stars_y[s])
    json.dumps(positions)
    
    #VISUALIZE
    #plt.close('all')
    #orbits
    colors = ['r', 'b', 'g', 'k', 'y', 'c', 'm',]
    plt.figure(figsize=(6,6))
    for s in range(number_of_stars):
        x = [positions[t][s][0] for t in range(timesteps)]
        y = [positions[t][s][1] for t in range(timesteps)]
        plt.scatter(x,y,c = colors[s%len(colors)], linewidth=0)
        plt.xlim((-width*1.1,width*1.1))
        plt.ylim((-width*1.1,width*1.1))
    #radius distribution
    #plt.figure()
    #plt.hist([stars_r[s] for s in range(number_of_stars)])
    #plt.xlabel('radius')
    #plt.ylabel('number of stars')

    
    