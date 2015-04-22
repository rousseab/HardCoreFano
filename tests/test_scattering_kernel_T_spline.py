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

eta = 1
hw_ext = 0.150

list_hw_ext = N.arange(0.100,0.220,0.001)

filename ='scattering_spline_TEST_T.nc'
n_xi_grid = 8001

SK = ScatteringKernel(mu, beta, kernel_Gamma_width, Green_Gamma_width, spline_order = 1)
SK.build_scattering_kernel_integrals(list_hw_ext, n_xi_grid = n_xi_grid )
SK.build_and_write_spline_parameters(filename)
SK.read_spline_parameters(filename)

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

list_xi_numeric = N.arange(-2,2,0.05)

list_re_TR_numeric, list_im_TR_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_TR, 'epsilon', global_parameters)

list_re_TR_der_numeric, list_im_TR_der_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_TR_derivative, 'epsilon', global_parameters)


list_re_TI_numeric, list_im_TI_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_TI, 'epsilon', global_parameters)

list_re_TI_der_numeric, list_im_TI_der_numeric  =\
    Drive_Numerical_Integrals(list_xi_numeric, integrand_TI_derivative, 'epsilon', global_parameters)


list_epsilon = N.arange(-1.5*D_cutoff,1.5*D_cutoff,0.010)

T_R = SK.get_spline_TR(list_epsilon, eta, hw_ext, sign_Gamma=1)
T_I = SK.get_spline_TI(list_epsilon, eta, hw_ext, sign_Gamma=1)



h = 8  # figure height
w = 18  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

list_ax = [ax1,ax2]

lw = 4

i_eta = SK.get_i_eta_index(eta)
i_hw_ext, i_hw_2, alpha_1,alpha_2 = SK.get_hw_indices_and_interpolation(hw_ext)
ax1.plot(SK.list_xi,N.real(SK.TR[i_eta,i_hw_ext,:]),'-',lw=lw,label='$T^R$')
ax1.plot(list_epsilon,N.real(T_R),'--',lw=lw,label='$T^R$ (spline)')
ax1.plot(list_xi_numeric ,list_re_TR_numeric,'o',label='NUMERICAL $T^R$')
ax1.plot(list_xi_numeric ,list_re_TR_der_numeric,'h--',label='Fermi derivative $T^R$')

ax1.plot(SK.list_xi,N.real(SK.TI[i_eta,i_hw_ext,:]),'-',lw=lw,label='$T^I$')
ax1.plot(list_epsilon,N.real(T_I),'--',lw=lw,label='$T^I$ (spline)')
ax1.plot(list_xi_numeric,list_re_TI_numeric,'o',label='NUMERICAL $T^I$')
ax1.plot(list_xi_numeric ,list_re_TI_der_numeric,'h--',label='NUMERICAL Fermi derivative $T^I$')


ax2.plot(SK.list_xi,N.imag(SK.TR[i_eta,i_hw_ext,:]),'-',lw=lw,label='$T^R$')
ax2.plot(list_epsilon,N.imag(T_R),'--',lw=lw,label='$T^R$ (spline)')
ax2.plot(list_xi_numeric ,list_im_TR_numeric,'o',label='NUMERICAL $T^R$')
ax2.plot(list_xi_numeric ,list_im_TR_der_numeric,'h--',label='NUMERICAL Fermi derivative $T^R$')


ax2.plot(SK.list_xi,N.imag(SK.TI[i_eta,i_hw_ext,:]),'-',lw=lw,label='$T^I$')
ax2.plot(list_epsilon,N.imag(T_I),'--',lw=lw,label='$T^I$ (spline)')
ax2.plot(list_xi_numeric,list_im_TI_numeric,'o',label='NUMERICAL $T^I$')
ax2.plot(list_xi_numeric,list_im_TI_der_numeric,'h--',label='NUMERICAL Fermi derivative $T^I$')



for ax in list_ax:
    ax.set_xlabel('$\epsilon$ (eV)')
    ax.set_xlim([-1,1])



ax1.set_ylabel('Re[$F(x)$]')
ax2.set_ylabel('Im[$F(x)$]')

# adjust the figure so that there is no wasted space


ax2.legend(  loc='center left', bbox_to_anchor=(1.00,0.50) , fancybox=True,shadow=True,  borderaxespad=0.)
fig1.suptitle('TESTING $T$: $\eta$ = %i, $\hbar\omega_{ext}$ = %i meV, $\mu$ = %i meV, $\Gamma_G$ = %i meV, $\Gamma_\gamma$ = %i meV'%(eta, 1000*hw_ext, 1000*mu,1000*Green_Gamma_width,1000*kernel_Gamma_width))

fig1.subplots_adjust(   right   =       0.70,
                        wspace  =       0.40,
                        hspace  =       0.20)



#plt.tight_layout()

plt.show()
