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

for delta in list_kernel_delta_wdith:


    list_z = list_xi+1j*delta

    """
    gamma  = get_gamma(list_z)
    Sigma_normalized = 0.5*gamma/D_cutoff*density
    re_Sigma = N.real(Sigma_normalized)
    im_Sigma = N.imag(Sigma_normalized)
    """

    normalization = 0.5/D_cutoff*density

    re_Sigma = get_KR(list_xi,delta)*normalization
    im_Sigma =-get_KI(list_xi,delta)*normalization



    line, = ax1.plot(list_xi/D_cutoff,-im_Sigma,'-',alpha=1.0,lw=lw,label='$\Gamma_\gamma= %4.1f$ meV'%(1000*delta))


    ax2.plot(list_xi/D_cutoff,re_Sigma,'-',alpha=1.0,lw=lw,label='$\Gamma_\gamma= %4.1f$ meV'%(1000*delta))



for ax in list_ax:
    ax.set_xlabel('$\\xi/D$')
    ax.legend(loc=0,fancybox=True,shadow=True)
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)

    ax.legend(loc=0,fancybox=True,shadow=True)


ax1.set_ylabel('-Im$[\Sigma(\\xi)]/(n_i D)$')
ax2.set_ylabel('Re$[\Sigma(\\xi)]/(n_i D)$')

# adjust the figure so that there is no wasted space
fig1.suptitle('Electronic Self-Energy, Self-Consistent Born approx. with broadening')
fig1.subplots_adjust(   left    =   0.07,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()

