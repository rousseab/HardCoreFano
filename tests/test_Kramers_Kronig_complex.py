#================================================================================
# Plot the Kramers-Kronig relations
#
#================================================================================

import common
reload(common)
from common import *

import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm


n_hw   = 2001
hw_max = 40.
hw_min =-40.

d_hw  = (hw_max-hw_min)/(n_hw-1.)

iloop  = N.arange(n_hw)
list_hw = hw_min+d_hw*iloop

Gamma = 0.010 

KKG = KramersKronig_Gamma(list_hw,Gamma)


fig = plt.figure(figsize=(16,10))

ax1  = fig.add_subplot(221)
ax2  = fig.add_subplot(222)
ax3  = fig.add_subplot(223)
ax4  = fig.add_subplot(224)

lw = 4

ax1.set_title('Localized kernel for KK integral')
ax2.set_title('Localized kernel for KK integral')
ax3.set_title('Real Part KK integral')
ax4.set_title('Imaginary Part KK integral')



ax1.set_ylabel('Re[F$(\omega)$]')
ax2.set_ylabel('Im[F$(\omega)$]')
ax3.set_ylabel('Re[KKF$(\omega)$]')
ax4.set_ylabel('Im[KKF$(\omega)$]')


delta = 1.0

# Lorentzian function
list_F = complex(1.,0.)/(list_hw+1j*delta)

#list_localized_kernel = 1j*N.imag(list_F)
list_localized_kernel = list_F

if Gamma > 0:
    list_exact_KK   = complex(2.,0.)/(list_hw+1j*Gamma+1j*delta)
elif Gamma < 0:
    list_exact_KK   = N.zeros_like(list_hw)

#list_Fr_exact = N.real(list_F_exact )

list_KK_F = KKG.apply_kramers_kronig_FFT_convolution(list_localized_kernel)

ax1.plot(list_hw,N.real(list_localized_kernel),'r-',lw=lw,label='Lorentzian')
ax2.plot(list_hw,N.imag(list_localized_kernel),'r-',lw=lw,label='Lorentzian')


ax3.plot(list_hw,N.real(list_KK_F ),'r-',lw=lw,label='Kramers-Kronig')
ax3.plot(list_hw,N.real(list_exact_KK),'g--',lw=lw,label='exact')


ax4.plot(list_hw,N.imag(list_KK_F),'r-',lw=lw,label='Kramers-Kronig')
ax4.plot(list_hw,N.imag(list_exact_KK),'g--',lw=lw,label='exact')


"""
ax2.plot(list_hw,N.real(list_Fr_w),'r-',lw=lw,label='Kramers-Kronig')
ax2.plot(list_hw,N.real(list_Fr_w_exact),'g--',lw=lw,label='exact')


fad = -1j*SS.wofz(list_hw)
list_Fi_w_exact  = 1j*N.imag(fad)
list_Fr_w_exact  =  N.real(fad)


#list_F_KK_w = KK.apply_kramers_kronig_FFT_convolution(list_Fi_w_exact)
list_Fi_w = KK.apply_kramers_kronig_FFT_convolution(list_Fi_w_exact)
list_Fr_w = KK.apply_kramers_kronig_FFT_convolution(list_Fi_w_exact)



ax3.plot(list_hw,N.imag(list_Fi_w),'r-',lw=lw,label='Kramers-Kronig')
ax3.plot(list_hw,N.imag(list_Fi_w_exact),'g--',lw=lw,label='exact')

"""


list_ax = [ax1,ax2,ax3,ax4]
for ax in list_ax:
        ax.grid(True,linestyle='-',color='grey',alpha=0.5)
        ax.set_xlabel('$\omega$')

fig.suptitle('Testing Kramers-Kronig: $\Gamma$ = %4.3f'%Gamma)

ax1.legend( loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)
ax3.legend( loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)

fig.subplots_adjust(    left    =       0.10,
                        bottom  =       0.10,
                        right   =       0.90,
                        top     =       0.90,
                        wspace  =       0.40,
                        hspace  =       0.40)



plt.show()
