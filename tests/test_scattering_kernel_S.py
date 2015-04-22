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

filename ='scattering_spline_TEST.nc'
n_xi_grid = 8001

SK = ScatteringKernel(mu, beta, kernel_Gamma_width, Green_Gamma_width, spline_order = 1)
SK.build_scattering_kernel_integrals(list_hw_ext, n_xi_grid = n_xi_grid )
#SK.read_spline_parameters(filename)

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


h = 8  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

list_ax = [ax1,ax2]

ax1.plot(SK.list_xi,N.real(SK.SR_inf ),'-',label='$S^\infty$')

ax1.plot(SK.list_xi,N.real(SK.SR-SK.SR_inf ),'-',label='$S^R$-$S^\infty$')

ax1.plot(list_xi_numeric ,list_re_SRinf_numeric,'o',label='NUMERICAL $S^R$-$S^\infty$')

ax1.plot(SK.list_xi,N.real(SK.SI),'-',label='$S^I$')
ax1.plot(list_xi_numeric,list_re_SI_numeric,'o',label='NUMERICAL $S^I$')


ax2.plot(SK.list_xi,N.imag(SK.SR_inf ),'-',label='$S^\infty$')
ax2.plot(SK.list_xi,N.imag(SK.SR-SK.SR_inf ),'-',label='$S^R$-$S^\infty$')

ax2.plot(list_xi_numeric ,list_im_SRinf_numeric,'o',label='NUMERICAL $S^R$-$S^\infty$')

ax2.plot(SK.list_xi,N.imag(SK.SI),'-',label='$S^I$')
ax2.plot(list_xi_numeric,list_im_SI_numeric,'o',label='NUMERICAL $S^I$')



for ax in list_ax:
    ax.set_xlabel('$x$ (eV)')


ax2.legend(loc=0)

ax1.set_ylabel('Re[$F(x)$]')
ax2.set_ylabel('Im[$F(x)$]')

# adjust the figure so that there is no wasted space


fig1.suptitle('TESTING $S$: $\mu$ = %i meV, $\Gamma_G$ = %i meV, $\Gamma_\gamma$ = %i meV'%(1000*mu,1000*Green_Gamma_width,1000*kernel_Gamma_width))
plt.tight_layout()

plt.show()
