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
from matplotlib import rcParams
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import json
import scipy.constants as cons

rcParams["savefig.dpi"] = 150
rcParams["axes.facecolor"] = '#000000'
rcParams["axes.linewidth"]=2
rcParams["axes.edgecolor"]='#dddddd'
rcParams["figure.facecolor"]   = '#000000'
rcParams["figure.edgecolor"]   = '#000000'
rcParams["axes.labelcolor"] = '#ffffff'
rcParams["xtick.color"] = '#ffffff'
rcParams["ytick.color"] = '#ffffff'


#define movement functions -- given an x,y return x',y'
def rotate(x, y, theta):
    """
    Rotate a cartesian coordinate given by *x*, *y* by an angle *theta*.
    Returns the new cartesian coordinate *x'*, *y'*.
    """
    return x*np.cos(theta) - y*np.sin(theta), y*np.cos(theta) + x*np.sin(theta)

def newton_angular_velocity(mass, radius):
    #print radius, np.sqrt(cons.G*mass/radius/ly2m), np.sqrt(cons.G*mass/radius/ly2m)/(ly2m*radius)
    return np.sqrt(cons.G*mass/radius/ly2m)
    
def iso_angular_velocity(Rc,Vinf, radius):
    tmp=1.0-Rc/radius/ly2m*np.arctan(radius*ly2m/Rc)
    #print radius, Vinf*np.sqrt(tmp), Vinf*np.sqrt(tmp) / (radius * ly2m)
    return Vinf*np.sqrt(tmp)
    
def nfw_angular_velocity(Vv,cv,Rv,radius):
    gcinv=np.log(1.0+cv)-1.0/(1.0/cv+1.0)
    s=radius/Rv/kpc2ly
    cvs=cv*s
    vsq=Vv*Vv/s/gcinv*(np.log(1.0+cvs)-1.0/(1.0/cvs+1))
    #print radius, np.sqrt(vsq), np.sqrt(vsq)/ (radius * ly2m)
    return np.sqrt(vsq)

#Will be called with input as JSON
if __name__ == "__main__":
    #set parameters
    #parameters = sys.argv[0]
    parameters = {"dark_matter":False, "model": "NFW", "amount_dark_matter": 6, "distribution": "Scenario A"}
    if not parameters['dark_matter']:
        pngnameroot= 'nodm'
    else:
        pngnameroot= parameters['model']        

    width = 50000
    number_of_stars = 200
    timesteps = 50
    lambda_ = 0.2  # exponential distribution of stars
    
    #constants
    dt = 1e14
    kpc2ly=3261.63344
    ly2m=1.0/1.05702341e-16
    H=2.1e-18
    kpc2m=3e19
    Msun=2e30 #kg
    gal='M33'
    Mv=2e11*Msun
    tmp=4.0/3.0*np.pi*2e-24
    Rv=np.power(Mv/tmp,0.33)*5
    M={'M33':1e40, 'Vinf':105400, 'Rc':1.39*kpc2m, 'Vv':H*10*Rv, 'cv':4.0,'Rv':Rv/kpc2m}
    
    
    #exponentially distribute stars from center of galaxy
    stars_r = [width*random.exponential(lambda_) for i in range(number_of_stars)]
    stars_theta = [random.uniform(2.0*np.pi) for i in range(number_of_stars)]
    stars_x = [stars_r[i]*np.cos(stars_theta[i]) for i in range(number_of_stars)]
    stars_y = [stars_r[i]*np.sin(stars_theta[i]) for i in range(number_of_stars)]
    

    
    #BUILD JSON
    positions = {0:{i:(stars_x[i], stars_y[i]) for i in range(number_of_stars)}} # save initial positions
    for t in range(1,timesteps):
        positions[t] = {}
        for s in range(number_of_stars):
            #rotate star's position
            if not parameters["dark_matter"]:
                #print "NEWTON"
                new_angle = dt*newton_angular_velocity(M[gal], stars_r[s]) / (stars_r[s] * ly2m)
            elif parameters["model"]=="ISO":
                #print "ISO"
                new_angle = dt*iso_angular_velocity(M['Rc'],M['Vinf'],stars_r[s]) / (stars_r[s] * ly2m)
            elif parameters["model"]=="NFW":
                #print "NFW"
                new_angle = dt* nfw_angular_velocity(M['Vv'],M['cv'],M['Rv'],stars_r[s]) / (stars_r[s] * ly2m)
            else:
                print "ERROR!"
            stars_x[s], stars_y[s] = rotate(stars_x[s], stars_y[s], new_angle)
            #save position
            positions[t][s] = (stars_x[s], stars_y[s])
    #json.dumps(positions)
    
    #VISUALIZE
    #plt.close('all')
    #orbits
    colors = ['r', 'b', 'g', 'k', 'c', 'm',]
    plt.figure(figsize=(6,6))
    for t in range(timesteps):
        for s in range(number_of_stars):
            x = positions[t][s][0]
            y = positions[t][s][1]

            plt.plot(x,y,'.',c = colors[s%len(colors)], alpha = t/float(timesteps), linewidth=0)
        plt.xlim(-30000,30000)
        plt.ylim(-30000,30000)
#        plt.savefig(pngnameroot+"%02d.png"%t)

    #
    fig = plt.figure(facecolor='k')
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-30000,30000), ylim=(-20000,30000))

    nodmplot, = ax.plot([], [],'.', marker='o',c = '#ffffff')#,label='NO DM')
    isoplot, = ax.plot([], [],'.', marker='o',c = '#ffffff')#,label='ISO')
    nfwplot, = ax.plot([], [],'.', marker='o',c = '#ffffff')#,label='NFW')

    def init():
        nodmplot.set_data([], [])
        isoplot.set_data([], [])
        nfwplot.set_data([], [])
        return nodmplot,isoplot,nfwplot
    
    def animate(t):
        x=[positions[t][s][0] for s in range(number_of_stars)]
        y=[positions[t][s][1] for s in range(number_of_stars)]
        if not parameters['dark_matter']:
            nodmplot.set_data(x,y)
        elif parameters['dark_matter'] and parameters['model']=='ISO':
            isoplot.set_data(x,y)
        elif parameters['dark_matter'] and parameters['model']=='NFW':
            nfwplot.set_data(x,y)
        return nodmplot,isoplot,nfwplot

    #add universe center
    plt.scatter([0], [0], c = 'y', marker='o', s=60)
    ani = animation.FuncAnimation(fig, animate, np.arange(1, timesteps),
                                 interval=25, blit=False, init_func=init)
    #set graph properties
    plt.xlim((-width*1.1,width*1.1))
    plt.ylim((-width*1.1,width*1.1))
    plt.show()
    
    
