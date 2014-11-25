import common
reload(common)
from common import *

# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm



kB = 8.6173324e-5 # eV/K
T  = 300 # K

beta = 1./(kB*T)

m  = N.arange(-50,50,1)
wm = N.pi/beta*(2*m+1)

mu = -0.4 # -400 meV
#mu = 0.4 # -400 meV
#mu = -0.0 # -400 meV

normalization = density/D_cutoff  

z     = 1j*wm+mu
gA    = get_gammaA(z)
gC    = get_gammaC(z)

xf    = N.pi/beta*(2*N.arange(-50.,50.,0.0051))
zf    = mu+1j*xf

gAf    = get_gammaA(zf)
gCf    = get_gammaC(zf)

delta = 0.050  # eV
gCf_real_axis = get_gammaC(xf+1j*delta)
gAf_real_axis = get_gammaA(xf+1j*delta)


h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(221)
ax2 = fig1.add_subplot(223)
ax3 = fig1.add_subplot(222)
ax4 = fig1.add_subplot(224)

list_ax = [ax1,ax2,ax3]

ms  = 10 # marker size
mew = 2  # marker edge width

y = N.imag(z)
yf= N.imag(zf)

for ax, ff in zip(list_ax,[gAf+gCf,gCf,gAf]):
    ax.plot(yf,normalization*N.real(ff),'r-',lw=1,label='__nolabel__')
    ax.plot(yf,normalization*N.imag(ff),'b-',lw=1,label='__nolabel__')
    

ax4.plot(xf,normalization*N.real(gAf_real_axis),'r--',lw=3,label='gA: Real Part')
ax4.plot(xf,normalization*N.imag(gAf_real_axis),'b--',lw=3,label='gA: Imaginary Part')
ax4.plot(xf,normalization*N.real(gCf_real_axis),'r-',lw=3,label='gC: Real Part')
ax4.plot(xf,normalization*N.imag(gCf_real_axis),'b-',lw=3,label='gC: Imaginary Part')




for ax, f in zip(list_ax,[gA+gC,gC,gA]):
    ax.plot(y,normalization*N.real(f),'ro',ms=8,mew=2,label='Real part')
    ax.plot(y,normalization*N.imag(f),'bo',ms=8,mew=2,label='Imaganary part')



for ax in list_ax:
    ax.set_xlabel('$\omega_m$ (eV)')
    ax.legend(loc=0,fancybox=True,shadow=True)
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)

ax4.legend(loc=0,fancybox=True,shadow=True)

ax4.set_xlabel('$\epsilon$ (eV)')
ax4.grid(True,linestyle='-',color='grey',alpha=0.5)

ax1.set_ylabel('$n \gamma(z)/D$')
ax2.set_ylabel('$n \gamma_C(z)/D$')
ax3.set_ylabel('$n \gamma_A(z)/D$')
ax4.set_ylabel('$n \gamma(x)/D$')

# adjust the figure so that there is no wasted space
fig1.suptitle('$\mu$ = %4.3f eV'%mu)
fig1.subplots_adjust(   left    =   0.07,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()
