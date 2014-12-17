import common
reload(common)
from common import *

import scipy.integrate as SI
# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

mpl.rcParams['font.size']='16'

def cutoff_denominator(list_x,tol):
    """
    This function returns an approximation to 1/(x+i delta) using the Faddeeva function
    """
    z   = list_x/tol

    den = -1j*SS.wofz(z)/tol*N.sqrt(N.pi)

    return den


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


    tol = 1e-6

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

mu    =-0.400  # -400 meV
#delta = 0.500 # 500 meV
kernel_Gamma_width = 0.200  # 200 meV

# fill the object with dummy variables
q_vector    = N.array([0.,0.])
list_hw     = N.array([0.])
delta_width = 0.

filename ='scattering_spline_mu=%4.3f_eV_delta=%4.3f_eV.nc'%(mu,kernel_Gamma_width)

os.system('cp %s %s'%(filename,'scattering_spline.nc'))

wedge = DummyWedge(N.array([q_vector]))

IG = IGridFunction( 'smooth', q_vector,list_hw, delta_width, kernel_Gamma_width, mu, beta, wedge)

normalization = density/D_cutoff  

filename ='scattering_spline_mu=%4.3f_eV_delta=%4.3f_eV.nc'%(mu,kernel_Gamma_width)

SK = ScatteringKernel(mu,beta,kernel_Gamma_width)
#SK.build_scattering_kernel()
#SK.build_and_write_spline_parameters(filename)
SK.read_spline_parameters(filename)

list_xi = N.arange(-2.*D_cutoff,2*D_cutoff,0.005)

h = 8   # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

list_ax = [ax1,ax2]

ms  = 10 # marker size
mew = 2  # marker edge width
lw  = 3


xi2=  -1.0

XI2 = [-1.,0.,1.]

hw = 0.200

for xi2 in XI2:
    list_xi2 = xi2*N.ones_like(list_xi)

    for eta,ax in zip([-1.,1.],list_ax):

        list_eta_hw = eta*N.array([hw])
        explicit_Y_smooth =  get_Y_smooth(list_xi, list_xi2, list_eta_hw, SK)[:,0]
        Y_smooth =  IG.get_Y_smooth(list_xi, list_xi2, list_eta_hw)[:,0]

        error = explicit_Y_smooth- Y_smooth 


        line, = ax.plot(list_xi,normalization*N.real(Y_smooth),'-',alpha=1.0,lw=lw,label='$\\xi_2 = %4.3f$ eV'%xi2)

        c=line.get_color()

        x = (xi2-eta*hw)*N.array([1,1])
        y = N.array(ax.get_ylim())
        ax.plot(x,y,c+'--',alpha=0.5,lw=lw,label='$\\xi = \\xi_2 - \eta \hbar\omega$')


ax1.set_title('$\eta =-1$')
ax2.set_title('$\eta =+1$')


for ax in list_ax:
    ax.set_xlabel('$\\xi$ (eV)')
    ax.legend(loc=0,fancybox=True,shadow=True)
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)

    ax.legend(loc=0,fancybox=True,shadow=True)


    ax.set_ylabel('$n/D Y(\\xi,\hbar\omega)$')

# adjust the figure so that there is no wasted space
fig1.suptitle('$\mu$ = %4.3f eV, $\hbar\omega$ = %4.3f eV, $\Gamma_\gamma$ = %4.3f eV'%(mu,hw,kernel_Gamma_width))
fig1.subplots_adjust(   left    =   0.07,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()

