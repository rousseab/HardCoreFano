#================================================================================
#
#
#   This script tests the behavior of the ill-behaved Kramers-Kronig integral
#   of f(xi)(K^R - K^inf), which only dies slowly with xi!
#
#================================================================================

import common
reload(common)
from common import *

# import the usual modules
from functools import partial
import scipy.integrate as SI
import numpy as N
import matplotlib.pyplot as plt
import seaborn; seaborn.set(font_scale=2)


# parameters
mu = -0.400 # eV
T = 300.
kernel_Gamma = 0.200 #eV 
Green_Gamma = 0.100 #eV 
beta = 1./(kB*T)

normalization = 2./Area/D_cutoff

h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(111)

fig2 = plt.figure(figsize=(w,h))
ax2 = fig2.add_subplot(221)
ax3 = fig2.add_subplot(222)
ax4 = fig2.add_subplot(223)
ax5 = fig2.add_subplot(224)


list_frac = [2.,4.,8]

n_xi_grid_0 = 1000 

def get_correction(list_xi,mu,Lambda_cutoff,Green_Gamma):

    
    prefactor = -1j*hvF**2/2.        
    den = 1./(list_xi+mu+1j*Green_Gamma)

    Term_1 = N.log( (Green_Gamma**2+(Lambda_cutoff+list_xi)**2)/(Lambda_cutoff-mu)**2)
    Term_2 = 1j*N.pi-1j*2*N.arctan((Lambda_cutoff+list_xi)/Green_Gamma)

    correction = prefactor*den*(Term_1+Term_2)

    return correction 

# Perform numerical integrals for a few values, to get an independent check

def get_integrand(xi, epsilon, mu, beta, Green_Gamma, kernel_Gamma):

    list_xi = N.array([xi])
    f = function_fermi_occupation(list_xi,0.,beta)[0]
    KR = get_KR(list_xi+mu,kernel_Gamma)[0]
    K_inf = get_KR_inf(list_xi+mu)[0]

    integrand = 1./(1j*N.pi)*f*(KR-K_inf)/(xi-epsilon-1j*Green_Gamma)

    return integrand 


list_xi_numeric = N.arange(-2.*D_cutoff,2.*D_cutoff,0.5)

parameters = {'mu':mu, 'beta':beta, 'Green_Gamma':Green_Gamma, 'kernel_Gamma':kernel_Gamma}
list_re_S_numeric, list_im_S_numeric =  Drive_Numerical_Integrals(list_xi_numeric, get_integrand,'epsilon', parameters)


for frac in list_frac:

    # build a regular linear energy grid
    n_xi_grid  =  frac*n_xi_grid_0 +1

    xi_max = frac*D_cutoff
    xi_min =-frac*D_cutoff
    d_xi   = (xi_max-xi_min)/(n_xi_grid-1.)
    iloop  = N.arange(n_xi_grid)

    list_xi = xi_min+d_xi*iloop

    # Compute the KK integral
    KK = KramersKronig_Gamma(list_xi,Green_Gamma)

    f = function_fermi_occupation(list_xi,0.,beta)
    KR = get_KR(list_xi+mu,kernel_Gamma)
    KR_inf = get_KR_inf(list_xi+mu)

    integral_kernel = f*(KR-KR_inf) 

    KK_SR_raw = KK.apply_kramers_kronig_FFT_convolution(integral_kernel)

    correction = get_correction(list_xi,mu,frac*D_cutoff,Green_Gamma)

    KK_SR = KK_SR_raw + correction

    ax1.plot(list_xi/D_cutoff,normalization*integral_kernel,label='max = %i D$_{cutoff}$'%frac)

    ax2.plot(list_xi/D_cutoff,normalization*N.real(KK_SR_raw),label='max = %i D$_{cutoff}$'%frac)
    ax3.plot(list_xi/D_cutoff,normalization*N.imag(KK_SR_raw),label='max = %i D$_{cutoff}$'%frac)
    ax4.plot(list_xi/D_cutoff,normalization*N.real(KK_SR),label='max = %i D$_{cutoff}$'%frac)
    ax5.plot(list_xi/D_cutoff,normalization*N.imag(KK_SR),label='max = %i D$_{cutoff}$'%frac)



ax4.plot(list_xi_numeric/D_cutoff,normalization*list_re_S_numeric,'o' ,label='NUMERIC')
ax5.plot(list_xi_numeric/D_cutoff,normalization*list_im_S_numeric,'o' ,label='NUMERIC')


for ax in [ax1,ax2,ax3,ax4,ax5]:
    ax.set_xlabel('$\\xi/D_{cutoff}$ ')

ax1.legend(loc=0,fancybox=True,shadow=True)
ax4.legend(loc=0,fancybox=True,shadow=True)

ax1.set_ylabel('$n/D$  ($K^R(\\xi)-K^{\infty}(\\xi)$)')

ax2.set_ylabel('$n/D$  Re[KK] (RAW)')
ax3.set_ylabel('$n/D$  Im[KK] (RAW)')
ax4.set_ylabel('$n/D$  Re[KK] (CORRECTED)')
ax5.set_ylabel('$n/D$  Im[KK] (CORRECTED)')



#ax1.set_ylim([-15.,15.])
#fig1.suptitle('Impurity Scattering kernel')

# adjust the figure so that there is no wasted space

plt.tight_layout()
plt.show()
