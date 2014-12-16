import common
reload(common)
from common import *

import os
#----------------------------------------
# import modules
#----------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties

import matplotlib.cm as cm
import matplotlib.tri as tri


def get_xi(k_vector,q_vector,mu,n1,n2,n3):

    eps_k  = N.real(function_epsilon_k(N.array([k_vector])))[0]
    eps_kq = N.real(function_epsilon_k(N.array([k_vector+q_vector])))[0]

    if n1 == 0:
        sign1 = -1.
    else:
        sign1 =  1.
    if n2 == 0:
        sign2 = -1.
    else:
        sign2 =  1.
    if n3 == 0:
        sign3 = -1.
    else:
        sign3 =  1.

    xi1 = sign1*eps_k-mu
    xi2 = sign2*eps_k-mu
    xi3 = sign3*eps_kq-mu

    return xi1, xi2, xi3


def build_matsubara_sum_I(k_vector,q_vector,mu,beta,n1,n2,n3,m_max, delta_width, list_z):
    xi1, xi2, xi3 = get_xi(k_vector,q_vector,mu,n1,n2,n3)

    # Define the Matsubara frequencies, symmetric about 0, all the way to the cutoff
    iwm = 1j*N.pi/beta*(2*N.arange(-m_max,m_max)+1) # in eV

    shift= 1j*delta_width*N.sign(-1j*iwm)

    iwm = iwm

    gamma = get_gamma(iwm+shift+mu) # in eV x a0^2

    G1 = complex(1.,0.)/(iwm-xi1) # in 1/eV
    G3 = complex(1.,0.)/(iwm-xi3) # in 1/eV

    gGG = gamma*G1*G3

    denominator = (iwm[:,N.newaxis]-xi2)**2-list_z[N.newaxis,:]**2
    numerator   = -2.*list_z[N.newaxis,:]

    M = numerator/denominator # in units of 1/eV

    # the summand is in units of  a0^2/eV
    summand = complex(1.,0.)/beta*M*gGG[:,N.newaxis]

    unit_conversion = Ha_to_eV
    I = N.sum(summand,axis=0)*unit_conversion  

    return I


def build_analytic_I_1(k_vector,q_vector,mu,beta,n1,n2,n3, delta_width, list_i_Omega_n):

    xi1, xi2, xi3 = get_xi(k_vector,q_vector,mu,n1,n2,n3)

    f1 = function_fermi_occupation(xi1,0,beta)
    f2 = function_fermi_occupation(xi2,0,beta)
    f3 = function_fermi_occupation(xi3,0,beta)

    g1 = N.real(get_gamma(xi1+mu+1j*delta_width)) # in eV x a0^2
    g3 = N.real(get_gamma(xi3+mu+1j*delta_width)) # in eV x a0^2

    L1 = f1*g1
    L3 = f3*g3

    I  = complex(0.,0.)*N.zeros_like(list_i_Omega_n) 

    for eta in [-1.,1.]:
    
        sign = -eta*N.sign(-1j*list_i_Omega_n)
        g2 = get_gamma(xi2+mu-eta*list_i_Omega_n+1j*delta_width*sign) # in eV x a0^2

        L2 = f2*g2

        Y12 = (L1-L2)/(xi1-xi2+eta*list_i_Omega_n)
        Y32 = (L3-L2)/(xi3-xi2+eta*list_i_Omega_n)

        I += eta*(Y12-Y32)/(xi1-xi3)


    unit_conversion = Ha_to_eV
    I = I*unit_conversion 

    return I


def build_analytic_I_2(k_vector,q_vector,mu,beta,n1,n2,n3, delta_width, list_i_Omega_n):


    xi1, xi2, xi3 = get_xi(k_vector,q_vector,mu,n1,n2,n3)

    list_xi, list_KK_ifKI = build_scattering_kernel(mu, beta, delta_width,abs_Omega_n=0.)

    L1 = interpolate_value(xi1, list_xi, list_KK_ifKI )
    L3 = interpolate_value(xi3, list_xi, list_KK_ifKI )

    I  = complex(0.,0.)*N.zeros_like(list_i_Omega_n) 

    for j, i_Omega_n in enumerate(list_i_Omega_n):

        abs_Omega_n = N.abs(i_Omega_n)
        sgn_Omega_n = N.sign(-1j*i_Omega_n)

        list_xi, list_KK_ifKI = build_scattering_kernel(mu, beta, delta_width,abs_Omega_n)

        L2_abs = interpolate_value(xi2, list_xi, list_KK_ifKI )
        L2_r = N.real(L2_abs) 
        L2_i = N.imag(L2_abs) 

        for eta in [-1.,1.]:
        
            L2 = complex(1.,0.)*L2_r-1j*eta*sgn_Omega_n*L2_i            

            Y12 = (L1-L2)/(xi1-xi2+eta*i_Omega_n)
            Y32 = (L3-L2)/(xi3-xi2+eta*i_Omega_n)

            I[j] += eta*(Y12-Y32)/(xi1-xi3)

    unit_conversion = Ha_to_eV
    I = I*unit_conversion 

    return I


def build_scattering_kernel(mu, beta, delta_width,abs_Omega_n):
    """
    Compute the Kramers-Kronig integral of the form   

    I = int_{-infty}^{infty} dx/(i pi) 1/(x-xi-i Omega_n) [ i f(x) K^I(x+mu) ]

    where K^I(x+mu ) = (gamma(x+mu+i delta) - gamma(x+mu-i delta))/(2i)
    """

    # Hardcode these parameters for now.

    # build a regular linear energy grid
    n_x   = 2001
    x_max = 2.*D_cutoff
    x_min =-2.*D_cutoff
    d_x   = (x_max-x_min)/(n_x-1.)
    iloop  = N.arange(n_x)
    list_x = x_min+d_x*iloop

    # Prepare the Kramers-Kronig object
    if abs_Omega_n < 1e-10:
        KK = KramersKronig(list_x)
    else:
        KK = KramersKronig_Gamma(list_x,abs_Omega_n)

    # Compute the anti-symmetric kernel
    g_plus    = get_gamma(list_x+mu+1j*delta_width)
    g_minus   = get_gamma(list_x+mu-1j*delta_width)

    KI = (g_minus-g_plus)/(2.*1j)

    # Compute the Fermi occupation factor. Note that the chemical potential should NOT be
    # passed to this function 
    f  = function_fermi_occupation(list_x,0.,beta)

    # the integrand in the Kramers-Kronig integral
    list_ifKI = 1j*f*KI

    # apply the Kramers-Kronig transformation to the anti-symmetric kernel
    list_KK_ifKI = KK.apply_kramers_kronig_FFT_convolution(list_ifKI)

    return list_x, list_KK_ifKI 


def interpolate_value(xi, list_xi, list_KK_ifKI ):

    # find the closest value in the array

    i_min = N.abs(list_xi-xi).argmin()

    if list_xi[i_min] >= xi:
        i1 = i_min-1
        i2 = i_min
    else:
        i1 = i_min
        i2 = i_min+1


    x1 = list_xi[i1]
    x2 = list_xi[i2]

    y1 = list_KK_ifKI[i1]
    y2 = list_KK_ifKI[i2]

    y  = (x2-xi)/(x2-x1)*y1+(xi-x1)/(x2-x1)*y2

    return y



#================================================================================
# Compute I
#================================================================================

T = 300. # K
beta = 1./(kB*T)
mu   = -0.400 # eV


k_vector = N.random.random(2)*twopia
q_vector = N.random.random(2)*twopia

#k_vector = N.array([ 0.63891072,  0.5939504 ])
#q_vector = N.array([ 0.11364837,  0.38050302])



wedge = DummyWedge(N.array([k_vector]))

n_max = 80

delta_width = 0.050
list_i_OMEGA_n = 1j*2*N.pi/beta*N.arange(-n_max,n_max+1)

n1 = 1
n2 = 1
n3 = 1

xi1, xi2, xi3 = get_xi(k_vector,q_vector,mu,n1,n2,n3)


m_max = 100000


fig = plt.figure(figsize=(12,8))
ax1  = fig.add_subplot(111)

list_ax = [ax1]
for ax in list_ax:
    ax.set_xlabel('$\Omega_n$ (eV)')

ax1.set_ylabel('Im[$I$] ($a_0^2/Ha$)')

list_delta_width = [0.050,0.150,0.200]
list_colors      = ['r','g','b']
list_symbols     = ['o','s','h']

for delta_width,c,s in zip(list_delta_width,list_colors,list_symbols): 
    explicit_I = build_matsubara_sum_I(k_vector,q_vector,mu,beta,n1,n2,n3,m_max, delta_width, list_i_OMEGA_n )

    ax1.plot(N.imag(list_i_OMEGA_n),N.imag(explicit_I),s+'-',c=c,ms=12,mew=2,label='$\delta$ = %4.3f eV'%delta_width)

    spline_file = 'scattering_spline_mu=%4.3f_eV_delta=%4.3f_eV.nc'%(mu,delta_width)
    os.system('cp %s scattering_spline.nc'%spline_file ) 

    analytic_I_1 = build_analytic_I_1(k_vector,q_vector,mu,beta,n1,n2,n3, delta_width, list_i_OMEGA_n )
    analytic_I_2 = build_analytic_I_2(k_vector,q_vector,mu,beta,n1,n2,n3, delta_width, list_i_OMEGA_n )

    analytic_I = analytic_I_1 + analytic_I_2 

    ax1.plot(N.imag(list_i_OMEGA_n),N.imag(analytic_I),s+'-',c=c,ms=4,mew=1,label='analytic')


ax1.grid(True,linestyle='-',color='grey',alpha=0.5)
ax1.legend(  loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)


"""
fig = plt.figure(figsize=(12,8))

ax1  = fig.add_subplot(131)
ax2  = fig.add_subplot(132)
ax3  = fig.add_subplot(133)

list_ax  = [ax1,ax2,ax3]

for ax in list_ax:
    ax.set_xlabel('$\Omega_n$ (eV)')

ax1.set_ylabel('Re[$I$] ($a_0^2/Ha$)')
ax2.set_ylabel('Im[$I$] ($a_0^2/Ha$)')


fig.suptitle('k = (%5.3f, %5.3f) $2\pi/a$, q = (%5.3f, %5.3f) $2\pi/a$'%(k_vector[0]/twopia, k_vector[1]/twopia, q_vector[0]/twopia, q_vector[1]/twopia))

ax1.plot(N.imag(list_i_OMEGA_n),N.real(explicit_I),'ro-',ms=8,mew=2,label='EXPLICIT')
ax2.plot(N.imag(list_i_OMEGA_n),N.imag(explicit_I),'ro-',ms=8,mew=2,label='EXPLICIT')


ax1.plot(N.imag(list_i_OMEGA_n),N.real(analytic_I),'ko-',ms=4,mew=2,label='analytic')
ax2.plot(N.imag(list_i_OMEGA_n),N.imag(analytic_I),'ko-',ms=4,mew=2,label='analytic')


ax3.plot(N.imag(list_i_OMEGA_n),N.imag(explicit_I-analytic_I),'ro-',ms=4,mew=2,label='difference')








ax1.legend(  loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)
ax2.legend(  loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)
ax3.legend(  loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)

for ax in [ax1,ax2,ax3]:
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)
"""


fig.subplots_adjust(    left    =       0.08,
                        bottom  =       0.16,
                        right   =       0.96,
                        top     =       0.90,
                        wspace  =       0.40,
                        hspace  =       0.40)


plt.show()
