import common
reload(common)
from common import *

import scipy.integrate as SI
import scipy.special  as  SS
# import the usual modules
import numpy as N
import matplotlib as mpl
import seaborn; seaborn.set(font_scale=2)
import matplotlib.pyplot as plt


kB = 8.6173324e-5 # eV/K
T  = 300 # K
beta = 1./(kB*T)

hw = 0.250 # eV

list_xi = N.arange(-1.,1.,0.001)

list_hw = [0.010,0.100,0.200]
list_c  = ['r','g','b']


h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(111)

list_ax = [ax1]

ms  = 10 # marker size
mew = 2  # marker edge width
lw  = 3

for hw,c in zip(list_hw,list_c):
    F2 = function_fermi_occupation(list_xi-hw,0.,beta)
    F1 = function_fermi_occupation(list_xi,0.,beta)

    dF = F2-F1
    ax1.plot(list_xi,dF,c+'-',lw=1.5*lw,label='$\hbar\omega$ = %3.1f meV'%(1000*hw),alpha=0.5)

    #dF_approx = -hw*d_Fermi_dxi(list_xi,0.,beta)+0.5*hw**2*d2_Fermi_dxi2(list_xi,0.,beta)
    dF_approx = -hw*d_Fermi_dxi(list_xi,0.,beta)

    ax1.plot(list_xi,dF_approx,c+'--',lw=lw,label='$\hbar\omega$ = %3.1f meV'%(1000*hw))

for ax in list_ax:
    ax.set_xlabel('$\\xi$ (eV)')
    ax.legend(loc=0,fancybox=True,shadow=True)

    ax.legend(loc=0,fancybox=True,shadow=True)


ax1.set_ylabel('$f(\\xi-\hbar\omega)-f(\\xi)$')

fig1.suptitle('Difference in Fermi-Dirac occupation factors')
# adjust the figure so that there is no wasted space

plt.show()
