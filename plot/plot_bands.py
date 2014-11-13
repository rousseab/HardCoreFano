import common
reload(common)
from common import *


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
# compute dispersion
#----------------------------------------

nkpoints  = 200
kpath     = Path(acell,nkpoints)
ksize     = kpath.list_k.shape[0]

epsilon_k = function_epsilon_k(kpath.list_k)


dist_max  = 0.8*N.linalg.norm(kpath.K-kpath.M)

dk_vec    = kpath.list_k-kpath.K
dk        = N.sqrt(dk_vec[:,0]**2+dk_vec[:,1]**2)
cone      = hvF*dk

restricted_cone = []
restricted_x    = []

for d, c,x in zip(dk, cone, kpath.list_x):
    if d < dist_max:
        restricted_cone.append(c)
        restricted_x.append(x)

restricted_cone = N.array(restricted_cone )
restricted_x    = N.array(restricted_x    )


#----------------------------------------
# plot dispersion
#----------------------------------------

fig = plt.figure(figsize=(14,10))
ax  = fig.add_subplot(111)


ax.plot(kpath.list_x, epsilon_k,'r-',lw=5)
ax.plot(kpath.list_x,-epsilon_k,'b-',lw=5)

ax.plot(restricted_x, restricted_cone,'g--',lw=5)
ax.plot(restricted_x,-restricted_cone,'g--',lw=5)



#ax.set_title(r'Tight-binding dispersion of Graphene')
ax.set_ylabel(r'$\epsilon_{\bf k}$ (eV)')

ax.set_xticks(kpath.list_xticks)
ax.set_xticklabels(kpath.list_labels)

ax.set_xlim(kpath.list_x[0],kpath.list_x[-1])

ax.grid(True,linestyle='-',color='grey',axis='y',alpha=0.5)
ax.grid(True,linestyle='-',color='black',axis='x',alpha=1.)

fig.subplots_adjust(    left    =       0.15,
                        bottom  =       0.10,
                        right   =       0.95,
                        top     =       0.92,
                        wspace  =       0.50,
                        hspace  =       0.50)

#plt.savefig('bands_graphene.pdf')
plt.show()
