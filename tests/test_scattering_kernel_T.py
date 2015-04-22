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
Green_Gamma_width  = 0.100   

list_hw_ext = N.arange(0.100,0.220,0.10)

SK = ScatteringKernel(mu, beta, kernel_Gamma_width, Green_Gamma_width, spline_order = 3)
SK.build_scattering_kernel_integrals(list_hw_ext, n_xi_grid = 4001)



eta    = 1 
i_eta  = N.argwhere(N.array(SK.list_eta) == eta)[0][0]

hw_ext = 0.200
i_hw_ext = N.argwhere(SK.list_hw_ext == hw_ext)[0][0]


N.where(eta == SK.list_eta)

# Compute the  kernels
def integrand_TR(xi, epsilon, eta, hw, mu, beta, kernel_Gamma_width, Green_Gamma_width):

    list_xi = N.array([xi])
    f     = function_fermi_occupation(list_xi,0.,beta)[0]
    f_ehw = function_fermi_occupation(list_xi+eta*hw,0.,beta)[0]
    KR = get_KR(list_xi+mu,kernel_Gamma_width)[0]

    den = 1./(xi-(epsilon+1j*Green_Gamma_width))

    integrand = 1./(1j*N.pi)*(f-f_ehw)*KR*den

    return integrand 

def integrand_TR_derivative(xi, epsilon, eta, hw, mu, beta, kernel_Gamma_width, Green_Gamma_width):

    list_xi = N.array([xi])

    df = d_Fermi_dxi(list_xi,0,beta)[0]

    KR = get_KR(list_xi+mu,kernel_Gamma_width)[0]

    den = 1./(xi-(epsilon+1j*Green_Gamma_width))

    integrand = 1./(1j*N.pi)*(-eta*hw*df)*KR*den

    return integrand 

def integrand_TI(xi, epsilon, eta, hw, mu, beta, kernel_Gamma_width, Green_Gamma_width):

    list_xi = N.array([xi])
    f     = function_fermi_occupation(list_xi,0.,beta)[0]
    f_ehw = function_fermi_occupation(list_xi+eta*hw,0.,beta)[0]
    KI = get_KI(list_xi+mu,kernel_Gamma_width)[0]

    den = 1./(xi-(epsilon+1j*Green_Gamma_width))

    integrand = 1./(1j*N.pi)*(f-f_ehw)*KI*den

    return integrand 

def integrand_TI_derivative(xi, epsilon, eta, hw, mu, beta, kernel_Gamma_width, Green_Gamma_width):

    list_xi = N.array([xi])
    KI = get_KI(list_xi+mu,kernel_Gamma_width)[0]

    den = 1./(xi-(epsilon+1j*Green_Gamma_width))

    df = d_Fermi_dxi(list_xi,0,beta)[0]

    integrand = -1./(1j*N.pi)*eta*hw*df*KI*den

    return integrand 


global_parameters = {'eta':eta, 'hw':hw_ext, 'mu':mu, 'beta':beta, 
                    'kernel_Gamma_width':kernel_Gamma_width, 
                    'Green_Gamma_width':Green_Gamma_width}

list_xi_numeric = N.arange(-2.,2.,0.25)

list_re_TR_numeric, list_im_TR_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_TR, 'epsilon', global_parameters)

list_re_TR_der_numeric, list_im_TR_der_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_TR_derivative, 'epsilon', global_parameters)


list_re_TI_numeric, list_im_TI_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_TI, 'epsilon', global_parameters)

list_re_TI_der_numeric, list_im_TI_der_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_TI_derivative, 'epsilon', global_parameters)


h = 8  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

list_ax = [ax1,ax2]

lw = 4
ms = 12

ax1.plot(SK.list_xi,N.real(SK.TR[i_eta,i_hw_ext,:]),'r-',lw=lw,label='$T^R$')
ax1.plot(list_xi_numeric ,list_re_TR_numeric,'o',ms=ms,label='NUMERICAL $T^R$')
ax1.plot(list_xi_numeric ,list_re_TR_der_numeric,'h--',ms=ms,label='Fermi derivative $T^R$')

ax1.plot(SK.list_xi,N.real(SK.TI[i_eta,i_hw_ext,:]),'b-',lw=lw,label='$T^I$')
ax1.plot(list_xi_numeric,list_re_TI_numeric,'o',ms=ms,label='NUMERICAL $T^I$')
ax1.plot(list_xi_numeric ,list_re_TI_der_numeric,'h--',ms=ms,label='NUMERICAL Fermi derivative $T^I$')


ax2.plot(SK.list_xi,N.imag(SK.TR[i_eta,i_hw_ext,:]),'r-',lw=lw,label='$T^R$')
ax2.plot(list_xi_numeric ,list_im_TR_numeric,'o',ms=ms,label='NUMERICAL $T^R$')
ax2.plot(list_xi_numeric ,list_im_TR_der_numeric,'h--',ms=ms,label='NUMERICAL Fermi derivative $T^R$')


ax2.plot(SK.list_xi,N.imag(SK.TI[i_eta,i_hw_ext,:]),'b-',lw=lw,label='$T^I$')
ax2.plot(list_xi_numeric,list_im_TI_numeric,'o',ms=ms,label='NUMERICAL $T^I$')
ax2.plot(list_xi_numeric,list_im_TI_der_numeric,'h--',ms=ms,label='NUMERICAL Fermi derivative $T^I$')





for ax in list_ax:
    ax.set_xlabel('$\epsilon$ (eV)')
    ax.set_xlim([-3,3])


ax2.legend(  loc='center left', bbox_to_anchor=(1.00,0.50) , fancybox=True,shadow=True,  borderaxespad=0.)

ax1.set_ylabel('Re[$F(x)$]')
ax2.set_ylabel('Im[$F(x)$]')


# adjust the figure so that there is no wasted space


fig1.suptitle('TESTING $T$: $\eta$ = %i, $\hbar\omega_{ext}$ = %i meV, $\mu$ = %i meV, $\Gamma_G$ = %i meV, $\Gamma_\gamma$ = %i meV'%(eta, 1000*hw_ext, 1000*mu,1000*Green_Gamma_width,1000*kernel_Gamma_width))
#plt.tight_layout()
fig1.subplots_adjust(   right =     0.70,
                        wspace=     0.40)


plt.show()
