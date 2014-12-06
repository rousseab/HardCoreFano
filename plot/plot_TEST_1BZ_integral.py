#================================================================================
# 
#       This script tests the convergence of the integral 
#
#           I1(mu) = 1/Omega sum_{k} delta(xi_k) 
#           I2(mu) = 1/Omega sum_{k} delta(xi_k) delta(xi_kq-hw) 
#
#       In the cone approximation, the result is given by 
#
#           I1_exact(mu) =  1/pi |mu|/(hvF)^2
#           I2_exact(mu) =  1/(2 pi^2) 1/(hvF)^2 (k_mu+k_w)/q 1/|sin(theta)|
#           
#                               cos( theta) = ( (k_mu+k_w)^2-k_mu^2-q^2)/(2 k_mu q)
#
#================================================================================
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

mu = 0.400 # eV
hw = 0.200

clip_energy = mu+0.5
grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma,\
                                clip_grid=True,clip_energy=clip_energy)

delta_width = 0.0250 # eV

K = 2./3.*twopia*N.array([1.,0.])

q  = N.array([0.200,0.])*twopia
nq = N.linalg.norm(q)
q_vec = q

I1_exact =  1./N.pi*N.abs(mu)/(hvF)**2.

k_mu = N.abs(mu)/hvF
k_w  = N.abs(hw)/hvF

cos_theta= ((k_mu+k_w)**2-k_mu**2-nq**2)/(2.*k_mu*nq)
theta = N.arccos(cos_theta)
I2_exact =  1./(2.*N.pi**2)/(hvF)**2*(k_mu+k_w)/nq/N.abs(N.sin(theta))

fig = plt.figure(figsize=(10,10))

ax  = fig.add_subplot(111)
ax.set_xlabel('$k_x$ ($2\pi/a$)')
ax.set_ylabel('$k_y$ ($2\pi/a$)')

def cutoff_denominator(list_x,fadeeva_width):
    """
    This function returns an approximation to 1/(x+i delta) using the Faddeeva function
    """
    z   = list_x/fadeeva_width

    den = -1j*SS.wofz(z)/fadeeva_width*N.sqrt(N.pi)

    return den

def delta_function(list_x,fadeeva_width):

    delta = -N.imag( cutoff_denominator(list_x,fadeeva_width) )/N.pi

    return delta

Integral1 = 0.
Integral2 = 0.

# 1/Omega sum_k = 1/(2pi)^2 int d2k

integration_factor = 1./(2.*N.pi)**2


list_cmax = []
list_cmin = []
list_images = []

for i,wedge in enumerate(grid.list_wedges):


    list_epsilon_k = function_epsilon_k(wedge.list_k)
    list_epsilon_kq= function_epsilon_k(wedge.list_k+q_vec)

    list_xik = list_epsilon_k-N.abs(mu) 
    list_xikq= list_epsilon_kq-N.abs(mu) 

    delta_xik = delta_function(list_xik,delta_width)
    delta_xikq= delta_function(list_xikq-hw,delta_width)

    list_Fk = delta_xik[:,N.newaxis]
    Integral1 += integration_factor*AreaIntegrator(wedge,list_Fk)[0]

    list_Fk = delta_xik[:,N.newaxis]*delta_xikq[:,N.newaxis]
    Integral2 += integration_factor*AreaIntegrator(wedge,list_Fk)[0]


    x = wedge.list_k[:,0]/twopia
    y = wedge.list_k[:,1]/twopia
    t = wedge.triangles_indices
    fc = (i+1)*N.arange(len(t))[::-1]
        
    #ax.triplot(x,y,triangles=t)
    #ax.tripcolor(x,y,triangles=t,facecolors=fc)
    image = ax.tripcolor(x,y, list_Fk[:,0], shading='gouraud', edgecolors='k',cmap=plt.cm.spectral_r)

    list_images.append(image)
    cmin, cmax = image.get_clim()

    list_cmax.append(cmax)
    list_cmin.append(cmin)

cmin = N.min(list_cmin)
cmax = N.max(list_cmax)

for image in list_images:
    image.set_clim([cmin,cmax])

axc = plt.colorbar(image)                                                                                           



# Plot Fermi surface

k_mu = N.abs(mu)/hvF
k_w  = N.abs(hw)/hvF

c = N.cos(N.pi/3.)
s = N.sin(N.pi/3.)
list_centers = N.array([ [-2./3.,0.], [ 2./3.,0.], 2./3.*N.array([c,s]), 
                           2./3.*N.array([c,-s]), 2./3.*N.array([-c,s]), 2./3.*N.array([-c,-s])]) 

rad_mu    = k_mu/twopia
rad_mu_w  = (k_mu+k_w)/twopia

th = N.arange(0.,2.*N.pi,0.01)
for center in list_centers:

    x = center[0]+rad_mu*N.cos(th)
    y = center[1]+rad_mu*N.sin(th)
    ax.plot(x, y,'k-',lw=2)


    x = center[0]-q[0]/twopia+rad_mu_w*N.cos(th)
    y = center[1]-q[1]/twopia+rad_mu_w*N.sin(th)
    ax.plot(x, y,'r-',lw=2)




ratio1 = N.real(Integral1)/I1_exact
ratio2 = N.real(Integral2)/I2_exact 

print 'ratio I1/I1_exact: %12.8e'%ratio1
print 'ratio I2/I2_exact: %12.8e'%ratio2


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
