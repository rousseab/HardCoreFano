#================================================================================
# 

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
import matplotlib.ticker as ticker

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

mpl.rcParams['font.size'] = 20.
legendfonts = FontProperties(size=16)

kB = 8.6173324e-5 # eV/K
T  = 300 # K
beta = 1./(kB*T)


nmax_coarse = 8
nmax_fine   = 256
n_blocks_coarse_to_fine = 2
#include_Gamma = True
include_Gamma = False

mu = -0.400 # eV
hw =  0.200

eta = -1.

t1 = time.time()
grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma,\
                                clip_grid=False,clip_energy=0.)
t2 = time.time()

print 'creating grid time: %5.3f sec.'%(t2-t1)


delta_width = 0.050 # eV

K = 2./3.*twopia*N.array([1.,0.])


fig1 = plt.figure(figsize=(14,10))

ax1 = fig1.add_subplot(111)
ax1.set_xlabel('$k_x$ ($2\pi/a$)')
ax1.set_ylabel('$k_y$ ($2\pi/a$)')


list_cmax_re = []
list_cmin_re = []
list_images_re = []

compute_time = 0.
for i,wedge in enumerate(grid.list_wedges):


    t1 = time.time()
    list_epsilon_k = function_epsilon_k(wedge.list_k)

    list_xi2k  = -list_epsilon_k-mu 


    Fermi2 = function_fermi_occupation(list_xi2k-eta*hw,0.,beta)
    Fermi1 = function_fermi_occupation(list_xi2k,0.,beta)

    integrand  = Fermi2-Fermi1


    #print 'max: %12.8e    %12.8e '%(units*N.abs(N.real(integrand)).max(),units*N.abs(N.imag(integrand)).max())

    list_Fk = integrand[:,N.newaxis]
    t2 = time.time()

    compute_time += t2-t1

    x = wedge.list_k[:,0]/twopia
    y = wedge.list_k[:,1]/twopia
    t = wedge.triangles_indices
    fc = (i+1)*N.arange(len(t))[::-1]
        
    image_re = ax1.tripcolor(x,y, N.real(list_Fk[:,0]), shading='gouraud', edgecolors='k',cmap=plt.cm.spectral_r)

    list_images_re.append(image_re)

    cmin_re, cmax_re = image_re.get_clim()

    list_cmax_re.append(cmax_re)
    list_cmin_re.append(cmin_re)



print 'computation time: %5.3f sec.'%compute_time


cmin_re = N.min(list_cmin_re)
cmax_re = N.max(list_cmax_re)


for image in list_images_re:
    image.set_clim([cmin_re,cmax_re])

axc_re = fig1.colorbar(image,format=ticker.FuncFormatter(fmt))


# Plot Fermi surface

k_w  = N.abs(hw)/hvF

c = N.cos(N.pi/3.)
s = N.sin(N.pi/3.)
list_centers = N.array([ [-2./3.,0.], [ 2./3.,0.], 2./3.*N.array([c,s]), 
                           2./3.*N.array([c,-s]), 2./3.*N.array([-c,s]), 2./3.*N.array([-c,-s])]) 

rad_mu    = N.abs(mu)/hvF/twopia

th = N.arange(0.,2.*N.pi,0.01)
for center in list_centers:

    x = center[0]+rad_mu*N.cos(th)
    y = center[1]+rad_mu*N.sin(th)

    ax1.plot(x, y,'k-',lw=2)




for line in ax1.xaxis.get_ticklines()+ax1.yaxis.get_ticklines():
    line.set_markeredgewidth(3)


ax1.set_xlim([-0.7,0.7])
ax1.set_ylim([-0.7,0.7])
fig1.gca().set_aspect('equal')

if eta == 1:
    title = '$f(\\xi-\hbar\omega) - f(\\xi)$'
else:
    title = '$f(\\xi+\hbar\omega) - f(\\xi)$'

fig1.suptitle(title+'\n$\mu$ = %3.1f meV, $\hbar\omega$ = %3.1f meV'%(1000*mu,1000*hw))

fig1.subplots_adjust(   left    =       0.10,
                        bottom  =       0.10,
                        right   =       0.90,
                        top     =       0.90,
                        wspace  =       0.20,
                        hspace  =       0.20)

plt.show()
