import common
reload(common)
from common import *

import scipy.special as SS

#----------------------------------------
# import modules
#----------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties

import matplotlib.cm as cm

mpl.rcParams['font.size'] = 50.
legendfonts = FontProperties(size=16)

#----------------------------------------
# function
#----------------------------------------

def cutoff_denominator(list_x):

    tol = 1e-3

    z   = list_x/tol

    #den = N.imag(SS.wofz(list_x))/tol
    den = N.imag(SS.wofz(z))/tol*N.sqrt(N.pi)

    return den


#----------------------------------------
# plot dispersion
#----------------------------------------

fig = plt.figure(figsize=(14,10))
ax1  = fig.add_subplot(211)
ax2  = fig.add_subplot(212)

list_x = N.arange(-5,5,0.0011)

f = cutoff_denominator(list_x)
g = 1./list_x

ax1.plot(list_x, f,'r-',lw=5)
ax1.plot(list_x, g,'b--',lw=5)
ax2.plot(list_x, list_x*f,'r-',lw=5)
ax2.plot(list_x, list_x*g,'b--',lw=5)

ax1.set_ylabel(r'1/x')
ax1.set_ylim([-4,4])

ax1.grid(True,linestyle='-',color='grey', alpha=0.5)
ax2.grid(True,linestyle='-',color='grey', alpha=0.5)

fig.subplots_adjust(    left    =       0.15,
                        bottom  =       0.10,
                        right   =       0.95,
                        top     =       0.92,
                        wspace  =       0.50,
                        hspace  =       0.50)

#plt.savefig('bands_graphene.pdf')
plt.show()
