import common
reload(common)
from common import *

import scipy.integrate as SI

from functools import partial

# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn;seaborn.set(font_scale=2)



kB = 8.6173324e-5 # eV/K
T  = 300 # K
beta = 1./(kB*T)

mu    = -0.400  
kernel_Gamma_width = 0.200  # 200 meV
Green_Gamma_width  = 0.050   

list_hw_ext = N.arange(0.100,0.220,0.001)
#list_hw_ext = N.arange(0.100,0.220,0.1)

filename ='scattering_spline_TEST.nc'
n_xi_grid = 8001

SK = ScatteringKernel(mu, beta, kernel_Gamma_width, Green_Gamma_width, spline_order = 1)
SK.build_scattering_kernel_integrals(list_hw_ext, n_xi_grid = n_xi_grid )
SK.read_spline_parameters(filename)

# Compute the  kernels
def integrand_SR_minus_S_infty(xi, epsilon, mu, beta, kernel_Gamma_width, Green_Gamma_width):

    list_xi = N.array([xi])
    f  = function_fermi_occupation(list_xi,0.,beta)[0]
    KR = get_KR(list_xi+mu,kernel_Gamma_width)[0]
    K_inf = get_KR_inf(list_xi+mu)

    integrand = 1./(1j*N.pi)*f*(KR-K_inf)/(xi-(epsilon+1j*Green_Gamma_width))

    return integrand 

def integrand_SI(xi, epsilon, mu, beta, kernel_Gamma_width, Green_Gamma_width):
    list_xi = N.array([xi])
    f  = function_fermi_occupation(list_xi,0.,beta)[0]
    KI = get_KI(list_xi+mu,kernel_Gamma_width)[0]

    integrand = 1./(1j*N.pi)*f*KI/(xi-(epsilon+1j*Green_Gamma_width))

    return integrand 


global_parameters = {'mu':mu, 'beta':beta, 
                    'kernel_Gamma_width':kernel_Gamma_width, 
                    'Green_Gamma_width':Green_Gamma_width}

list_xi_numeric = N.arange(-D_cutoff,D_cutoff,0.5)


list_re_SRinf_numeric, list_im_SRinf_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_SR_minus_S_infty, 'epsilon', global_parameters)

list_re_SI_numeric, list_im_SI_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_SI, 'epsilon', global_parameters)

SR_inf = get_SR_inf(list_xi_numeric,mu,Green_Gamma_width)


list_epsilon = N.arange(-1.5*D_cutoff,1.5*D_cutoff,0.010)

S_R = SK.get_spline_SR(list_epsilon, sign_Gamma=1)
S_I = SK.get_spline_SI(list_epsilon, sign_Gamma=1)

h = 8  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

list_ax = [ax1,ax2]

lw=4

ax1.plot(SK.list_xi,N.real(SK.SR),'-',lw=lw,label='$S^R$')
ax1.plot(list_epsilon,N.real(S_R),'--',lw=lw,label='$S^R$ (spline)')
ax1.plot(list_xi_numeric ,list_re_SRinf_numeric+N.real(SR_inf),'o',label='NUMERICAL $S^R$')

ax1.plot(SK.list_xi,N.real(SK.SI),'-',lw=lw,label='$S^I$')
ax1.plot(list_epsilon,N.real(S_I),'--',lw=lw,label='$S^I$ (spline)')
ax1.plot(list_xi_numeric,list_re_SI_numeric,'o',label='NUMERICAL $S^I$')

ax2.plot(SK.list_xi,N.imag(SK.SR),'-',lw=lw,label='$S^R$')
ax2.plot(list_epsilon,N.imag(S_R),'--',lw=lw,label='$S^R$ (spline)')
ax2.plot(list_xi_numeric ,list_im_SRinf_numeric+N.imag(SR_inf),'o',label='NUMERICAL $S^R$')

ax2.plot(SK.list_xi,N.imag(SK.SI),'-',lw=lw,label='$S^I$')
ax2.plot(list_epsilon,N.imag(S_I),'--',lw=lw,label='$S^I$ (spline)')
ax2.plot(list_xi_numeric,list_im_SI_numeric,'o',label='NUMERICAL $S^I$')



for ax in list_ax:
    ax.set_xlabel('$x$ (eV)')


ax2.legend(loc=0)

ax1.set_ylabel('Re[$F(x)$]')
ax2.set_ylabel('Im[$F(x)$]')

# adjust the figure so that there is no wasted space


fig1.suptitle('TESTING $S$: $\mu$ = %i meV, $\Gamma_G$ = %i meV, $\Gamma_\gamma$ = %i meV'%(1000*mu,1000*Green_Gamma_width,1000*kernel_Gamma_width))
#plt.tight_layout()

ax2.legend(  loc='center left', bbox_to_anchor=(1.00,0.50) , fancybox=True,shadow=True,  borderaxespad=0.)

fig1.subplots_adjust(   right   =       0.70,
                        wspace  =       0.40,
                        hspace  =       0.20)

plt.show()
