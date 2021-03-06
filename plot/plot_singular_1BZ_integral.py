#================================================================================
# 
#       This script tests the convergence of the singular term,
#
#       I_singular = -sum_eta eta [ f(xi2-eta hw) - f(xi2)] KR(xi2+mu-eta hw) 
#                                   1/(xi1-xi2+eta hw + i eta delta) 
#                                   1/(xi3-xi2+eta hw + i eta delta) 
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

clip_energy = N.abs(mu)+0.5

t1 = time.time()
grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma,\
                                clip_grid=True,clip_energy=clip_energy)
t2 = time.time()

print 'creating grid time: %5.3f sec.'%(t2-t1)


delta_width = 0.050 # eV


c  = N.cos(N.pi/3.)
s  = N.sin(N.pi/3.)

K = 2./3.*twopia*N.array([c,s])
q  = K-0.5*N.array([c,s])*twopia

nq = N.linalg.norm(q)
q_vec = q
hatQ  = q/nq
hatzQ = N.cross(N.array([0.,0.,1.]),hatQ)[:2]

fig1 = plt.figure(figsize=(14,10))
fig2 = plt.figure(figsize=(14,10))

ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
list_ax = [ax1,ax2]
for ax in list_ax:
    ax.set_xlabel('$k_x$ ($2\pi/a$)')
    ax.set_ylabel('$k_y$ ($2\pi/a$)')

def cutoff_denominator(list_x,fadeeva_width):
    """
    This function returns an approximation to 1/(x+i delta) using the Faddeeva function
    """
    z   = list_x/fadeeva_width

    den = -1j*SS.wofz(z)/fadeeva_width*N.sqrt(N.pi)

    return den


# 1/Omega sum_k = 1/(2pi)^2 int d2k
integration_factor = 1./(2.*N.pi)**2

# fundamental units conversion
units = Ha_to_eV

list_cmax_re = []
list_cmin_re = []

list_cmax_im = []
list_cmin_im = []

list_images_re = []
list_images_im = []

Integral = complex(0.,0.)

sign_1 = -1.
sign_2 = -1.
sign_3 = -1.

compute_time = 0.
for i,wedge in enumerate(grid.list_wedges):


    t1 = time.time()
    list_epsilon_k = function_epsilon_k(wedge.list_k)
    list_epsilon_kq= function_epsilon_k(wedge.list_k+q_vec)

    list_xi1k  = sign_1*list_epsilon_k-mu 
    list_xi2k  = sign_2*list_epsilon_k-mu 
    list_xi3kq = sign_3*list_epsilon_kq-mu


    integrand = complex(0.,0.)*N.zeros_like(list_xi1k)
    for eta in [-1.,1.]:

        KR = get_KR(list_xi2k+mu-eta*hw,delta_width)

        Fermi2 = function_fermi_occupation(list_xi2k-eta*hw,0.,beta)
        Fermi1 = function_fermi_occupation(list_xi2k,0.,beta)

        den12  = cutoff_denominator(list_xi1k  - list_xi2k+eta*hw,eta*delta_width)
        den32  = cutoff_denominator(list_xi3kq - list_xi2k+eta*hw,eta*delta_width)

        integrand  += -eta*(Fermi2-Fermi1)*KR*den12*den32


    #print 'max: %12.8e    %12.8e '%(units*N.abs(N.real(integrand)).max(),units*N.abs(N.imag(integrand)).max())

    list_Fk = units*integrand[:,N.newaxis]
    Integral += integration_factor*AreaIntegrator(wedge,list_Fk)[0]
    t2 = time.time()
    compute_time += t2-t1

    x = wedge.list_k[:,0]/twopia
    y = wedge.list_k[:,1]/twopia
    t = wedge.triangles_indices
    fc = (i+1)*N.arange(len(t))[::-1]
        
    image_re = ax1.tripcolor(x,y, N.real(list_Fk[:,0]), shading='gouraud', edgecolors='k',cmap=plt.cm.spectral_r)
    image_im = ax2.tripcolor(x,y, N.imag(list_Fk[:,0]), shading='gouraud', edgecolors='k',cmap=plt.cm.spectral_r)

    list_images_re.append(image_re)
    list_images_im.append(image_im)

    cmin_re, cmax_re = image_re.get_clim()
    cmin_im, cmax_im = image_im.get_clim()

    list_cmax_re.append(cmax_re)
    list_cmin_re.append(cmin_re)

    list_cmax_im.append(cmax_im)
    list_cmin_im.append(cmin_im)


print 'computation time: %5.3f sec.'%compute_time


cmin_re = N.min(list_cmin_re)
cmax_re = N.max(list_cmax_re)

cmin_im = N.min(list_cmin_im)
cmax_im = N.max(list_cmax_im)

for image in list_images_re:
    image.set_clim([cmin_re,cmax_re])

axc_re = fig1.colorbar(image,format=ticker.FuncFormatter(fmt))


for image in list_images_im:
    image.set_clim([cmin_im,cmax_im])

axc_im = fig2.colorbar(image,format=ticker.FuncFormatter(fmt))


axc_re.set_label(r'$a_0^2/Ha$')
axc_im.set_label(r'$a_0^2/Ha$')

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

    for ax in list_ax:
        ax.plot(x, y,'k-',lw=2)


    for eta in [-1.,1]:
        c = N.cos(th)
        s = N.sin(th)
        K = ((hw/hvF)**2-nq**2)/2./(nq*c-eta*hw/hvF)

        indk = N.where( (K >= 0)*(hvF*K < 1.0) )[0]

        K = K[indk]
        c = c[indk]
        s = s[indk]
        check = N.sqrt(K**2+nq**2+2.*K*nq*c)-K-eta*hw/hvF

        x = center[0]+K/twopia*(hatQ[0]*c+hatzQ[0]*s)
        y = center[1]+K/twopia*(hatQ[1]*c+hatzQ[1]*s)
        for ax in list_ax:
            ax.plot(x, y,'r--',lw=2)


print '  I          : %12.8e   %+12.8e j'%(N.real(Integral),N.imag(Integral))


for ax in list_ax:
    for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
        line.set_markeredgewidth(3)


    ax.set_xlim([-0.7,0.7])
    ax.set_ylim([-0.7,0.7])
    fig1.gca().set_aspect('equal')
    fig2.gca().set_aspect('equal')


title = 'Singular $I$, q = (%3.2f, %3.2f) 2 $\pi/a$'%(q[0]/twopia,q[1]/twopia)

fig1.suptitle(title+'\n Real Part')
fig2.suptitle(title+'\n Imaginary Part')

for fig in [fig1,fig2]:
    fig.subplots_adjust(    left    =       0.10,
                            bottom  =       0.10,
                            right   =       0.90,
                            top     =       0.90,
                            wspace  =       0.20,
                            hspace  =       0.20)

plt.show()
