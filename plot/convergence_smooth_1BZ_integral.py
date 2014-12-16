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

def cutoff_denominator(list_x,fadeeva_width):
    """
    This function returns an approximation to 1/(x+i delta) using the Faddeeva function
    """
    z   = list_x/fadeeva_width

    den = -1j*SS.wofz(z)/fadeeva_width*N.sqrt(N.pi)

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

q  = N.array([0.005,0.005])*twopia
nq = N.linalg.norm(q)
q_vec = q

include_Gamma = False

list_nmax_coarse = [16,16,16,32,32,32,64,64,64,128,128,128]
list_nmax_fine   = [64,128,256,64,128,256,64,128,256,128,256,512]

list_n_blocks_coarse_to_fine = [4,4,4,8,8,8,16,16,16,32,32,32]

mu = -0.400 # eV
hw =  0.200
delta_width = 0.050 # eV



filename = 'scattering_spline_mu=%4.3f_eV_delta=%4.3f_eV.nc'%(mu,delta_width)
SK = ScatteringKernel(mu,beta,delta_width)
#SK.build_scattering_kernel()
#SK.build_and_write_spline_parameters(filename)
SK.read_spline_parameters(filename)


# 1/Omega sum_k = 1/(2pi)^2 int d2k
integration_factor = 1./(2.*N.pi)**2

# fundamental units conversion
units = Ha_to_eV

sign_1 = -1.
sign_2 = -1.
sign_3 = -1.


print '#                                            Convergence parameters '
print '#       I                     coarse  fine  blocks   k-points  grid time (sec)   compute time (sec)'
for nmax_coarse, nmax_fine, n_blocks_coarse_to_fine in zip(list_nmax_coarse,list_nmax_fine,list_n_blocks_coarse_to_fine ): 



    t1 = time.time()
    grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma, clip_grid=False)
    t2 = time.time()


    grid_time = t2-t1

    nkpoints = 12*grid.list_wedges[0].list_k.shape[0]

    Integral = complex(0.,0.)
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



        list_Fk = units*integrand[:,N.newaxis]
        Integral += integration_factor*AreaIntegrator(wedge,list_Fk)[0]
        t2 = time.time()
        compute_time += t2-t1



    print '  %12.8e               %3i    %3i      %2i     %6i      %5.3f            %5.3f'%(N.real(Integral), nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,nkpoints,grid_time,compute_time)
