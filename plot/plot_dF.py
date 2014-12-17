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
import matplotlib.tri as tri

#================================================================================
# Compute dF
#================================================================================

T = 300. # K
beta = 1./(kB*T)
mu   = -0.400 # eV

nmax_coarse = 24
nmax_fine   = 120
n_blocks_coarse_to_fine = 5

include_Gamma = False
grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma )

hw = 0.150 # eV
k_mu     = N.abs(mu)/hvF
k_w      = N.abs(hw)/hvF

K = twopia*2./3*N.array([1.,0.])

q_vector = k_mu*N.array([0.5,0.00])


fig = plt.figure(figsize=(12,10))

ax  = fig.add_subplot(111)


ax.set_xlabel('$k_x$ ($2\pi/a$)')

mymap = plt.get_cmap("spectral")

list_min = []
list_max = []

list_images = []

for i,wedge in enumerate(grid.list_wedges):
        
    dF = dF_GridFunction(q_vector, mu, beta, wedge)

    S  = dF.df

    x = wedge.list_k[:,0]/twopia
    y = wedge.list_k[:,1]/twopia

    t = wedge.triangles_indices

    fc = S[t].mean(axis=1)
    image = ax.tripcolor(x,y,triangles=t,facecolors=fc,cmap=plt.cm.spectral_r,  edgecolors='w')
    list_images.append(image)

    list_min.append(N.min(fc))
    list_max.append(N.max(fc))

vmin = N.min(list_min)
vmax = N.max(list_max)

norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
for image in list_images:
        image.set_norm(norm)



fig.colorbar(image,ax=ax)


c = N.cos(N.pi/3.)
s = N.sin(N.pi/3.)
list_centers = N.array([ [0.,0.], [-2./3.,0.], [ 2./3.,0.], 2./3.*N.array([c,s]),
                       2./3.*N.array([c,-s]), 2./3.*N.array([-c,s]), 2./3.*N.array([-c,-s])])

rad_max = (2*k_mu+k_w)/twopia
rad_mid = (2*k_mu-k_w)/twopia
rad_min = k_w/twopia

rad_F = k_mu/twopia
list_centers_F = N.array([ [-2./3.,0.], [ 2./3.,0.], 2./3.*N.array([c,s]),
                       2./3.*N.array([c,-s]), 2./3.*N.array([-c,s]), 2./3.*N.array([-c,-s])])


th = N.arange(0.,2.*N.pi,0.01)

for center in list_centers_F:

    x = center[0]+rad_F*N.cos(th)
    y = center[1]+rad_F*N.sin(th)

    ax.plot(x, y,'g-',lw=2)


for center in list_centers:

    x = center[0]+rad_max*N.cos(th)
    y = center[1]+rad_max*N.sin(th)

    ax.plot(x, y,'k-',lw=2)

    x = center[0]+rad_min*N.cos(th)
    y = center[1]+rad_min*N.sin(th)

    ax.plot(x, y,'k-',lw=2)


for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
    line.set_markeredgewidth(3)


ax.set_xlim([-0.7,0.7])
ax.set_ylim([-0.7,0.7])

fig.suptitle('dF')


fig.gca().set_aspect('equal')



fig.subplots_adjust(    left    =       0.16,
                        bottom  =       0.16,
                        right   =       0.96,
                        top     =       0.98,
                        wspace  =       0.20,
                        hspace  =       0.20)


plt.show()
