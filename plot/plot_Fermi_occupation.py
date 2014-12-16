import common
reload(common)
from common import *

import scipy.integrate as SI
import scipy.special  as  SS
# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm


kB = 8.6173324e-5 # eV/K
T  = 300 # K
beta = 1./(kB*T)

hw = 0.250 # eV

list_xi = N.arange(-1.,1.,0.001)

list_hw = [-0.300, 0.010,0.200,0.300]


h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(111)

list_ax = [ax1]

ms  = 10 # marker size
mew = 2  # marker edge width
lw  = 3

for hw in list_hw:
    F2 = function_fermi_occupation(list_xi-hw,0.,beta)
    F1 = function_fermi_occupation(list_xi,0.,beta)

    dF = F2-F1
    ax1.plot(list_xi,dF,'-',lw=lw,label='$\hbar\omega$ = %3.1f meV'%(1000*hw))



for ax in list_ax:
    ax.set_xlabel('$\\xi$ (eV)')
    ax.legend(loc=0,fancybox=True,shadow=True)
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)

    ax.legend(loc=0,fancybox=True,shadow=True)


ax1.set_ylabel('$f(\\xi-\hbar\omega)-f(\\xi)$')

fig1.suptitle('Difference in Fermi-Dirac occupation factors')
# adjust the figure so that there is no wasted space
fig1.subplots_adjust(   left    =   0.12,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()
