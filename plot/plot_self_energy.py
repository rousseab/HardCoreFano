#================================================================================
#
# This script plots the self energy 
#
#================================================================================

import common
reload(common)
from common import *

import scipy.integrate as SI
# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

mpl.rcParams['font.size']='20'


list_xi = N.arange(-5.0*D_cutoff,5.0*D_cutoff,0.005)

h = 12   # figure height
w = 20  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

list_ax = [ax1,ax2]

ms  = 10 # marker size
mew = 2  # marker edge width
lw  = 3



list_kernel_delta_wdith = [0.05,0.100,0.150,0.200]

normalization = 0.5/D_cutoff*density
for delta in list_kernel_delta_wdith:

    re_Sigma = get_KR(list_xi,delta)*normalization
    im_Sigma =-get_KI(list_xi,delta)*normalization

    line, = ax1.plot(list_xi/D_cutoff,-im_Sigma,'-',alpha=1.0,lw=lw,label='$\Gamma_\gamma= %4.1f$ meV'%(1000*delta))
    ax2.plot(list_xi/D_cutoff,re_Sigma,'-',alpha=1.0,lw=lw,label='__nolabel__')


I0 = 10. # eV a0^2

for I0 in [1.,10.,50.]:
    re_Sigma_I0 = I0*normalization
    ax2.plot(list_xi/D_cutoff,re_Sigma_I0*N.ones_like(list_xi),'--',alpha=1.0,lw=lw,label='$I_0= %4.1f$ eV $a_0^2$'%(I0))




for ax in list_ax:
    ax.set_xlabel('$\\xi/D$',labelpad=20)
    ax.legend(loc=0,fancybox=True,shadow=True)
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)

    ax.legend(loc=0,fancybox=True,shadow=True)


ax1.set_ylabel(r'-$N/N_{imp}$ $\times$ Im$[\Sigma(\xi)]/D$')
ax2.set_ylabel(r'$N/N_{imp}$ $\times$ Re$[\Sigma(\xi)]/D$')

ax1.set_xlim([-0.8,0.8])
ax2.set_xlim([-0.8,0.8])

# adjust the figure so that there is no wasted space
#fig1.suptitle('Electronic Self-Energy, Self-Consistent Born approx. with broadening')
fig1.subplots_adjust(   left    =   0.07,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()

