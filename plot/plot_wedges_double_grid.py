#================================================================================
# Plot a nice tesselation of the 1BZ
#
#================================================================================
import common
reload(common)
from common import *

import time
#----------------------------------------
# import modules
#----------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties

import matplotlib.cm as cm

mpl.rcParams['font.size'] = 20.
legendfonts = FontProperties(size=16)

nmax_coarse = 32
nmax_fine   = 192
n_blocks_coarse_to_fine = 5

nmax_coarse = 8
nmax_fine   = 512
n_blocks_coarse_to_fine = 2

#include_Gamma = True
include_Gamma = False

mu = -0.400 # eV
hw =  0.250 # eV


clip_energy = N.abs(mu)+0.5 # eV

#grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma )
t1=time.time()
grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma , \
                                clip_grid=True,clip_energy=clip_energy)
t2=time.time()
print 'time: %12.4f'%(t2-t1)


NK = 12*len(grid.list_wedges[0].list_k)

print 'number of k vectors: %i'%NK

dE_c = N.sqrt(grid.fx**2+grid.fy**2)/nmax_coarse*hvF
dE_f = N.sqrt(grid.fx**2+grid.fy**2)/nmax_fine*hvF

tf = grid.irreducible_wedge.triangles_indices[-1]
print '# Energy resolution for the grid '
print '#    sub-grid        hvF dk  (meV)'
print '#---------------------------------'

print '     coarse           %4.1f    '%(1000*dE_c)
print '      fine            %4.1f    '%(1000*dE_f)


fig = plt.figure(figsize=(10,10))

ax  = fig.add_subplot(111)
ax.set_xlabel('$k_x$ ($2\pi/a$)')
ax.set_ylabel('$k_y$ ($2\pi/a$)')


Integral = 0.


for i,wedge in enumerate(grid.list_wedges):

    """
    # try to clip
    fy = N.sqrt(3.)/3*twopia # 
    fx = fy/N.sqrt(3.)
    K_point = N.array([fx,fy]) 

    energies = hvF*N.sqrt(N.sum((wedge.list_k-K_point)**2,axis=1))

    map = N.nan*N.ones(len(energies))

    I = N.where( energies <  clip_energy )[0]
    map[I] = N.arange(len(I))

    new_list_k = wedge.list_k[I]

    truth_table = False*N.ones(len(energies),dtype=bool)
    truth_table[I] = True 

    truth_triangle = truth_table[wedge.triangles_indices]
    truth_triangle = truth_triangle[:,0]*truth_triangle[:,1]*truth_triangle[:,2]

    It = N.where(truth_triangle)[0]
    new_Jacobians      = wedge.Jacobians[It]
    new_cross_products = wedge.cross_products[It]

    new_triangle_indices = N.array(map[wedge.triangles_indices[It]],dtype=int)
    """

    


    ones = N.ones(len(wedge.list_k))
    list_Fk = ones[:,N.newaxis]
    Integral += AreaIntegrator(wedge,list_Fk)[0]


    x = wedge.list_k[:,0]/twopia
    y = wedge.list_k[:,1]/twopia
    t = wedge.triangles_indices
    fc = (i+1)*N.arange(len(t))[::-1]
        

    ax.triplot(x,y,triangles=t)
    #ax.tripcolor(x,y,triangles=t,facecolors=fc)


# Plot Fermi surface

k_mu = N.abs(mu)/hvF
k_hw = N.abs(hw)/hvF

c = N.cos(N.pi/3.)
s = N.sin(N.pi/3.)
list_centers = N.array([ [-2./3.,0.], [ 2./3.,0.], 2./3.*N.array([c,s]), 
                           2./3.*N.array([c,-s]), 2./3.*N.array([-c,s]), 2./3.*N.array([-c,-s])]) 

rad_mu  = k_mu/twopia

rad_mu_plus_hw  = (k_mu+k_hw)/twopia
rad_mu_minus_hw = (k_mu-k_hw)/twopia


rad_upper_limit = clip_energy/hvF/twopia

th = N.arange(0.,2.*N.pi,0.01)
for center in list_centers:

    x = center[0]+rad_mu*N.cos(th)
    y = center[1]+rad_mu*N.sin(th)
    ax.plot(x, y,'k-',lw=2)

    x = center[0]+rad_upper_limit*N.cos(th)
    y = center[1]+rad_upper_limit*N.sin(th)
    ax.plot(x, y,'g--',lw=2)




    for rad in [rad_mu_plus_hw ,rad_mu_minus_hw ]: 
        x = center[0]+rad*N.cos(th)
        y = center[1]+rad*N.sin(th)
        ax.plot(x, y,'r-',lw=2)



FBZ_area = (2.*N.pi)**2/Area
error = N.abs(FBZ_area - Integral)

print 'Integration error : %12.8e'%error

for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
    line.set_markeredgewidth(3)


ax.set_xlim([-0.7,0.7])
ax.set_ylim([-0.7,0.7])
fig.gca().set_aspect('equal')
fig.suptitle('First Brillouin Zone Tesselation')


fig.subplots_adjust(    left    =       0.20,
                        bottom  =       0.20,
                        right   =       0.90,
                        top     =       0.90,
                        wspace  =       0.20,
                        hspace  =       0.20)

plt.show()
