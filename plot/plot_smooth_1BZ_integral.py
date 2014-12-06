#================================================================================
# 
#       This script tests the convergence of the smooth term,
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



mpl.rcParams['font.size'] = 20.
legendfonts = FontProperties(size=16)


def get_Y_smooth(list_xi1, list_xi2, list_eta_hw, SK):
    """
    This routine computes the "smooth" part of the Y function,

    Y(xi1,xi2,eta_hw) =  ( L(xi1,0) -L(xi2-eta_hw, 0 ) ) / (xi1-[xi2-eta_hw]).

    It is understood that as xi1 -> xi2-eta_hw, the function should remain smooth. 
    Care must be taken to avoid singularities.

    We'll assume that list_xi1 and list_xi2 have the same dimension [nk], and list_eta_hw is [nw]
    """

    L1  = SK.get_L_oneD(list_xi1)[:,N.newaxis]

    u   = list_xi2[:,N.newaxis]-list_eta_hw[N.newaxis,:]

    shape = u.shape

    L2  =  SK.get_L_oneD(u.flatten()).reshape(shape)


    numerator = L1-L2

    denominator = list_xi1[:,N.newaxis]-list_xi2[:,N.newaxis]+list_eta_hw[N.newaxis,:]


    # Where the denominator is too close to zero, replace by value for argument close to zero.


    tol = 1e-5

    I,J = N.nonzero( N.abs(denominator) < tol)
    
    # substitute with safe value, to not divide by zero        
    denominator[I,J] = 1.
    Y_smooth = numerator/denominator

    # replace by derivative 

    L1_safe_plus  = SK.get_L_oneD(list_xi1[I]+tol)
    L1_safe_minus = SK.get_L_oneD(list_xi1[I]-tol)

    # note: u[I,J] is a 1D array, not a 2D array. This is because I and J are arrays, not slicing operators
    
    L2_safe_plus  = SK.get_L_oneD(u[I,J]+tol)
    L2_safe_minus = SK.get_L_oneD(u[I,J]-tol)


    # I'm not sure why this strange derivative works better, but I find that it does for tol large (for debugging)
    Y_smooth[I,J] = 0.5*((L1_safe_plus+L2_safe_plus)-(L1_safe_minus+L2_safe_minus))/(2.*tol)


    return Y_smooth


kB = 8.6173324e-5 # eV/K
T  = 300 # K
beta = 1./(kB*T)


nmax_coarse = 32
#nmax_fine   = 256
nmax_fine   = 128
n_blocks_coarse_to_fine = 8
#include_Gamma = True
include_Gamma = False

mu = -0.400 # eV
hw =  0.200
delta_width = 0.050 # eV

t1 = time.time()
grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma, clip_grid=False)
t2 = time.time()

print 'creating grid time: %5.3f sec.'%(t2-t1)



filename = 'scattering_spline_mu=%4.3f_eV_delta=%4.3f_eV.nc'%(mu,delta_width)
SK = ScatteringKernel(mu,beta,delta_width)
#SK.build_scattering_kernel()
#SK.build_and_write_spline_parameters(filename)
SK.read_spline_parameters(filename)



K = 2./3.*twopia*N.array([1.,0.])

q  = N.array([0.0250,0.045])*twopia
nq = N.linalg.norm(q)
q_vec = q

fig1 = plt.figure(figsize=(10,10))

ax1 = fig1.add_subplot(111)
list_ax = [ax1]
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


list_images_re = []

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

        list_eta_hw = N.array([eta*hw])

        Y12 = get_Y_smooth(list_xi1k, list_xi2k, list_eta_hw, SK)[:,0]
        Y32 = get_Y_smooth(list_xi3kq, list_xi2k, list_eta_hw, SK)[:,0]


        integrand  += eta*(Y12-Y32)/(list_xi1k-list_xi3kq)


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

    list_images_re.append(image_re)

    cmin_re, cmax_re = image_re.get_clim()

    list_cmax_re.append(cmax_re)
    list_cmin_re.append(cmin_re)



print 'computation time: %5.3f sec.'%compute_time


cmin_re = N.min(list_cmin_re)
cmax_re = N.max(list_cmax_re)


for image in list_images_re:
    image.set_clim([cmin_re,cmax_re])

axc_re = fig1.colorbar(image)


# Plot Fermi surface


c = N.cos(N.pi/3.)
s = N.sin(N.pi/3.)
list_centers = N.array([ [-2./3.,0.], [ 2./3.,0.], 2./3.*N.array([c,s]), 
                           2./3.*N.array([c,-s]), 2./3.*N.array([-c,s]), 2./3.*N.array([-c,-s])]) 

rad_mu    = N.abs(mu)/hvF/twopia
rad_hw    = N.abs(hw)/hvF/twopia

th = N.arange(0.,2.*N.pi,0.01)
for center in list_centers:

    x = center[0]+rad_mu*N.cos(th)
    y = center[1]+rad_mu*N.sin(th)

    for ax in list_ax:
        ax.plot(x, y,'k-',lw=2)


    for eta in [-1,1]:
        c = N.cos(th)
        s = N.sin(th)

        x = center[0]+(rad_mu+eta*rad_hw)*N.cos(th)
        y = center[1]+(rad_mu+eta*rad_hw)*N.sin(th)
        for ax in list_ax:
            ax.plot(x, y,'r-',lw=2)


print '  I          : %12.8e   %+12.8e j'%(N.real(Integral),N.imag(Integral))


for ax in list_ax:
    for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
        line.set_markeredgewidth(3)


    ax.set_xlim([-0.7,0.7])
    ax.set_ylim([-0.7,0.7])
    fig1.gca().set_aspect('equal')

fig1.suptitle('Real Integrand')

fig1.subplots_adjust(   left    =       0.20,
                        bottom  =       0.20,
                        right   =       0.90,
                        top     =       0.90,
                        wspace  =       0.20,
                        hspace  =       0.20)

plt.show()
