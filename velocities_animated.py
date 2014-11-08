import os,sys
import numpy as np
import pylab as pl
import scipy.constants as cons
#%pylab inline
from matplotlib import rcParams
import matplotlib.animation as animation
rcParams["savefig.dpi"] = 150
rcParams["savefig.dpi"] = 150
rcParams["axes.facecolor"] = '#000000'
rcParams["axes.linewidth"]=5
rcParams["axes.edgecolor"]='#dddddd'
rcParams["figure.facecolor"]   = '#000000'
rcParams["figure.edgecolor"]   = '#000000'
rcParams["axes.labelcolor"] = '#ffffff'
rcParams["xtick.color"] = '#ffffff'
rcParams["ytick.color"] = '#ffffff'
kpc2ly=3261.63344
ly2m=1.0/1.05702341e-16
H=2.1e-18
kpc2m=3e19
Msun=2e30 #kg

def newton(M,R):
    return np.sqrt(cons.G*M/R/ly2m)
def iso(Rc,Vinf,R):
    tmp=1.0-Rc/R/ly2m*np.arctan(R*ly2m/Rc)
    return Vinf*np.sqrt(tmp)
def nfw(Vv,cv,Rv,R):
    gcinv=np.log(1.0+cv)-1.0/(1.0/cv+1.0)
    s=R/Rv/kpc2ly
    cvs=cv*s
    vsq=Vv*Vv/s/gcinv*(np.log(1.0+cvs)-1.0/(1.0/cvs+1))
    return np.sqrt(vsq)

#finite radii
gal='M33'
Mv=2e11*Msun
tmp=4.0/3.0*np.pi*2e-24
Rv=np.power(Mv/tmp,0.33)*5
M={'M33':1e40, 'Vinf':105400, 'Rc':1.39*kpc2m, 'Vv':H*10*Rv, 'cv':4.0,'Rv':Rv/kpc2m}
R=np.arange(1,40)*1000.
model={}



model['nodm']={'v':newton(M[gal],R)}

model['ISO']={'v':iso(M['Rc'],M['Vinf'],R)}

model['NFW']={'v':nfw(M['Vv'],M['cv'],M['Rv'],R)}
      


fig = pl.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(4000,39000), ylim=(0,150000))

nodmplot, = ax.plot([], [],label='NO DM')
isoplot, = ax.plot([], [],label='ISO')
nfwplot, = ax.plot([], [],label='NFW')
#pl.plot(R,model['nodm']['v'],label='NO DM')
#pl.plot(R,model['ISO']['v'],label='ISO')
#pl.plot(R,model['NFW']['v'],label='NFW')
l=pl.legend()

dt = 0.05
t = np.arange(0.0, 39/dt, dt)

def init():
    nodmplot.set_data([], [])
    isoplot.set_data([], [])
    nfwplot.set_data([], [])
    return nodmplot,isoplot,nfwplot

def animate(i):
    nodmplot.set_data(R[:i],model['nodm']['v'][:i])
    isoplot.set_data(R[:i],model['ISO']['v'][:i])
    nfwplot.set_data(R[:i],model['NFW']['v'][:i])
    return nodmplot,isoplot,nfwplot
ani = animation.FuncAnimation(fig, animate, np.arange(1, len(t)),
    interval=25, blit=False, init_func=init)

for text in l.get_texts():
    text.set_color("white")
pl.xlim(4000,39000)
pl.xlabel("distance to galaxy center (ly)")
pl.ylabel("rotational velocity (km/s)")
pl.show()
