#================================================================================
#
# This script plots the self energy 
#
#================================================================================

import common
reload(common)
from common import *

# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from matplotlib.font_manager import  FontProperties
import matplotlib.cm as cm
legendfonts = FontProperties(size=16)
from mpltools import color


mpl.rcParams['font.size']= 20
#mpl.rcParams['text.usetex'] = True



h = 8   # figure height
w = 8  # figure width

fig1 = plt.figure(figsize=(w,h))
ax = fig1.add_subplot(111)

lw  = 6

hw = N.arange(-0.1,2.,0.001) # eV
hw_meV = 1000*hw

T = 300. # K

beta = 1./(kB*T)

mu = -0.45 # eV


sigma_interband = 1./8.*( N.tanh(0.25*beta*(hw+2*mu)) +  N.tanh(0.25*beta*(hw-2*mu)))
ax.plot(hw_meV,sigma_interband ,'--',c='b',alpha=1.0,lw=lw,label='INTERBAND $\mu = %i$ meV'%(1000*mu),zorder=5)


for Gamma in [0.025,0.100,0.200,0.250]:


    # Form in Heinz paper
    # sigma_intra = 2.*N.log(2)/N.pi/beta*Gamma/(hw**2+Gamma**2)

    # Form in  Carbotte's paper, Eq. 14 (our Gamma is 2xTheir Gamma) 
    sigma_intra = 2.*N.log(2*N.cosh(0.5*beta*mu))/N.pi/beta*Gamma/(hw**2+Gamma**2)

    ax.plot(hw_meV,sigma_interband+sigma_intra,'-',alpha=1.0,lw=lw,label='$\Gamma = %i$ meV'%(1000*Gamma))



ax.set_ylabel(r'$\sigma(\omega)$ ($e^2/\hbar$)')

ax.set_xlabel(r'$\hbar\omega$ (meV)')

#ax.set_xlim([0.,1500])
ax.set_xlim([0.,900])
ax.set_ylim([0.,0.8])

ax.legend(loc=1,fancybox=True,shadow=True,prop=legendfonts)

ax.tick_params('both', length=12, width=2, which='major',direction='out')
ax.tick_params('both',  length=8, width=1, which='minor',direction='out')

# adjust the figure so that there is no wasted space
fig1.subplots_adjust(   left    =   0.15,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   0.30, 
                        hspace  =   None)

plt.show()

