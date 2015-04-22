import common
reload(common)
from common import *


# import the usual modules
import numpy as N
import matplotlib.pyplot as plt
import seaborn; seaborn.set(font_scale=2)


mu = -0.4
list_Gamma_kernel = [0.025,0.05,0.150,0.200]

frac = 2.0
list_xi = N.arange(-frac*D_cutoff,frac*D_cutoff,0.01)

h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)


KR_linear = -2*N.pi*(hvF)**2/D_cutoff*list_xi/D_cutoff

for Gamma_kernel,ls in zip(list_Gamma_kernel,['-','-','--','-']):

    KR = get_KR(list_xi,Gamma_kernel)
    KR_inf = get_KR_inf(list_xi)

    normalization = 2./Area/D_cutoff


    ax1.plot(list_xi,normalization*KR ,label='$\Gamma_\gamma$ = %4.1f meV'%(1000*Gamma_kernel))
    ax2.plot(list_xi,normalization*(KR-KR_linear) ,label='$\Gamma_\gamma$ = %4.1f meV'%(1000*Gamma_kernel))

ax1.plot(list_xi,normalization*KR_linear,'k-',label='Linear approximation')
ax1.plot([mu,mu],[-6,6],'k--',label='$\mu$= %3.1f meV'%(1000*mu))

ax1.set_xlim([-frac*D_cutoff,frac*D_cutoff])

ax1.set_xlabel('$\\xi$ (eV)')
ax1.set_ylabel('$n/D$  $K^R(\\xi)$')
ax1.legend(loc=0,fancybox=True,shadow=True)



ax1.set_ylim([-15.,15.])
fig1.suptitle('Impurity Scattering kernel')

# adjust the figure so that there is no wasted space

plt.tight_layout()
plt.show()
