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


n_hw   = 1001
hw_max = 10.
hw_min =-10.

d_hw  = (hw_max-hw_min)/(n_hw-1.)

iloop  = N.arange(n_hw)
list_hw = hw_min+d_hw*iloop


Gamma = 0.010 

KK = KramersKronig(list_hw)
#KK = KramersKronig_Gamma(list_hw,Gamma)

#KKP = KramersKronig_KP(list_hw,Gamma)
#KKd = KramersKronig_Kd(list_hw,Gamma)


fig = plt.figure(figsize=(16,10))

ax1  = fig.add_subplot(221)
ax2  = fig.add_subplot(222)
ax3  = fig.add_subplot(223)
ax4  = fig.add_subplot(224)

lw = 4

ax1.set_title('Imaginary Part, Lorentzian')
ax2.set_title('Real Part, Lorentzian')
ax3.set_title('Imaginary Part, Faddeeva')
ax4.set_title('Real Part, Faddeeva')



ax1.set_ylabel('Im[f$(\omega)$]')

ax2.set_ylabel('Re[f$(\omega)$]')
ax3.set_ylabel('Im[f$(\omega)$]')
ax4.set_ylabel('Re[f$(\omega)$]')



delta = 1.0

list_Fi_w_exact = -complex(0.,1.)*delta/(delta**2+list_hw**2)
list_Fr_w_exact =  complex(1.,0.)*list_hw/(delta**2+list_hw**2)

list_Fr_w   = KK.apply_kramers_kronig_FFT_convolution(list_Fi_w_exact)
list_Fi_w   = KK.apply_kramers_kronig_FFT_convolution(list_Fr_w_exact)

#list_Fi_w = KKd.apply_kramers_kronig_FFT_convolution(list_Fi_w_exact)
#list_Fr_w = KKP.apply_kramers_kronig_FFT_convolution(list_Fi_w_exact)

ax1.plot(list_hw,N.imag(list_Fi_w),'r-',lw=lw,label='Kramers-Kronig')
ax1.plot(list_hw,N.imag(list_Fi_w_exact),'g--',lw=lw,label='exact')

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

ax4.plot(list_hw,N.real(list_Fr_w),'r-',lw=lw,label='Kramers-Kronig')
ax4.plot(list_hw,N.real(list_Fr_w_exact),'g--',lw=lw,label='exact')


list_ax = [ax1,ax2,ax3,ax4]
for ax in list_ax:
        ax.grid(True,linestyle='-',color='grey',alpha=0.5)
        ax.set_xlabel('$\omega$')

fig.suptitle('Testing Kramers-Kronig: $\Gamma$ = %4.3f'%Gamma)

ax1.legend( loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)
ax2.legend( loc = 0, fancybox=True,shadow=True,  borderaxespad=0.)

fig.subplots_adjust(    left    =       0.10,
                        bottom  =       0.10,
                        right   =       0.90,
                        top     =       0.90,
                        wspace  =       0.40,
                        hspace  =       0.40)



plt.show()
