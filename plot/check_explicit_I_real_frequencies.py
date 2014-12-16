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



def build_matsubara_sum_I(k_vector,q_vector,mu,beta,n1,n2,n3,m_max, delta_width, list_z):

        eps_k  = function_epsilon_k(N.array([k_vector]))[0]
        eps_kq = function_epsilon_k(N.array([k_vector+q_vector]))[0]

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

        # Define the Matsubara frequencies, symmetric about 0, all the way to the cutoff
        iwm = 1j*N.pi/beta*(2*N.arange(-m_max,m_max)+1) # in eV

        shift = 1j*delta_width*N.sign(-1j*iwm)
        gamma = get_gamma(iwm+mu+shift) # in eV x a0^2

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

        return I, xi1,xi2,xi3




#================================================================================
# Compute I
#================================================================================

#k_vector = N.random.random(2)*twopia
#q_vector = N.random.random(2)*twopia

#k_vector = N.array([ 0.63891072,  0.5939504 ])
#q_vector = N.array([ 0.11364837,  0.38050302])



list_hw  = N.arange(0.001,15.0,0.010) # eV


n1 = 0
n2 = 0
n3 = 0



T = 300. # K
beta = 1./(kB*T)
mu   = -0.400 # eV

k_vector = 1.00001*N.abs(mu)/hvF*N.array([1.,0.])
q_vector = N.abs(mu)/hvF*N.array([-2.,0.])

wedge = DummyWedge(N.array([k_vector]))
#q_vector = N.random.random(2)*twopia


m_max = 10000

delta_width = 0.050 # eV


spline_file = 'scattering_spline_mu=%4.3f_eV_delta=%4.3f_eV.nc'%(mu,delta_width)
os.system('cp %s scattering_spline.nc'%spline_file ) 


list_z = list_hw #+1j*delta_width
explicit_I,xi1,xi2,xi3 = build_matsubara_sum_I(k_vector,q_vector,mu,beta,n1,n2,n3,m_max, delta_width, list_z)


I_smooth   = IGridFunction('smooth', q_vector,list_hw, delta_width, mu, beta, wedge)
I_singular = IGridFunction('singular', q_vector,list_hw, delta_width, mu, beta, wedge)

key   = (n1,n2,n3)
index = I_smooth.index_dictionary[key]

Ism = I_smooth.I[index][0,:]
Isi = I_singular.I[index][0,:]
I   = Ism+Isi


fig = plt.figure(figsize=(12,8))

ax1  = fig.add_subplot(121)
ax2  = fig.add_subplot(122)

list_ax  = [ax1,ax2]

for ax in list_ax:
    ax.set_xlabel('$\hbar\omega$ (eV)')

ax1.set_ylabel('Re[$I$] ($a_0^2/Ha$)')
ax2.set_ylabel('Im[$I$] ($a_0^2/Ha$)')


fig.suptitle('k = (%5.3f, %5.3f) $2\pi/a$, q = (%5.3f, %5.3f) $2\pi/a$'%(k_vector[0]/twopia, k_vector[1]/twopia, q_vector[0]/twopia, q_vector[1]/twopia))

ax1.plot(list_hw,N.real(Ism),'r-',lw=4,label='smooth')
ax2.plot(list_hw,N.imag(Ism),'r-',lw=4,label='smooth')

ax1.plot(list_hw,N.real(Isi),'g-',lw=4,label='singular')
ax2.plot(list_hw,N.imag(Isi),'g-',lw=4,label='singular')

ax1.plot(list_hw,N.real(I),'b-',lw=4,label='sum')
ax2.plot(list_hw,N.imag(I),'b-',lw=4,label='sum')





ax1.plot(list_hw,N.real(explicit_I),'k--',lw=4,label='EXPLICIT')
ax2.plot(list_hw,N.imag(explicit_I),'k--',lw=4,label='EXPLICIT')


ax1.legend(  loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)
ax2.legend(  loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)

fig.subplots_adjust(    left    =       0.08,
                        bottom  =       0.16,
                        right   =       0.96,
                        top     =       0.90,
                        wspace  =       0.40,
                        hspace  =       0.40)


plt.show()
