import common
reload(common)
from common import *


# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm


mu = -0.4
list_delta = [0.025,0.05,0.150,0.200]

frac = 0.4
list_xi = N.arange(-frac*D_cutoff,frac*D_cutoff,0.01)

h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(111)



ms  = 10 # marker size
mew = 2  # marker edge width


ms  = 8
mew = 2
lw  = 4

for delta,ls in zip(list_delta,['-','-','--','-']):

    KR = get_KR(list_xi,delta)

    normalization = 2./Area/D_cutoff


    ax1.plot(list_xi,normalization*KR,ls=ls,lw=lw,label='$\delta$ = %4.1f meV'%(1000*delta))



KR_linear = -2*N.pi*(hvF)**2/D_cutoff*list_xi/D_cutoff

ax1.plot(list_xi,normalization*KR_linear,'k-',ls=ls,lw=lw,label='Linear approximation')

ax1.plot([mu,mu],[-6,6],'k--',lw=lw,label='$\mu$= %3.1f meV'%(1000*mu))

ax1.set_xlim([-frac*D_cutoff,frac*D_cutoff])

ax1.set_xlabel('$\\xi$ (eV)')
ax1.set_ylabel('$n/D$  $K^R(\\xi)$')
ax1.legend(loc=0,fancybox=True,shadow=True)
ax1.grid(True,linestyle='-',color='grey',alpha=0.5)



ax1.set_ylim([-15.,15.])
fig1.suptitle('Impurity Scattering kernel')

# adjust the figure so that there is no wasted space
plt.subplots_adjust(    left    =   0.10,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()
