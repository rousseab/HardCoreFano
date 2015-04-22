#================================================================================
# 
#       This script tests AreaIntegrator method.
#
#
#================================================================================
import common
reload(common)
from common import *

import scipy.integrate as SI

def identity_function(list_k):
    """
    Return identity
    """
    return complex(1.,0.)*N.ones(len(list_k))
def get_exact_identity_result():

    First_Brillouin_Zone_Area = (2.*N.pi)**2/Area

    return First_Brillouin_Zone_Area 


def gaussian_function(list_k):
    """
    Return Gaussian which dies off quickly at the edge.
    """
    x = list_k[:,0]
    y = list_k[:,1]

    sig_x  = 0.25
    sig_y  = 0.15

    exponent = -0.5*x**2/sig_x**2-0.5*y**2/sig_y**2
    den   = 1./(2*N.pi*sig_x*sig_y)

    return den*N.exp(exponent)

def get_exact_gaussian_result():
    """
    Gaussian chosen to integrate to 1.
    """
    test_G = lambda x,y: gaussian_function(N.array([[x,y]]))
    gfun = lambda x: -N.infty
    hfun = lambda x: N.infty
    A = SI.dblquad(test_G, -N.infty, N.infty, gfun,hfun)

    return A[0]


nmax_coarse = 64
nmax_fine   = 256
n_blocks_coarse_to_fine = 24
include_Gamma = True

grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma, clip_grid=False)

list_hw = complex(1.,0.)*N.random.random(10)

print '#=============================================================='
print '# Function       exact result     Area integrator    diff' 
print '#=============================================================='

list_functions = [identity_function,gaussian_function]
list_exact     = [get_exact_identity_result(),get_exact_gaussian_result()]
list_names     = ['constant','Normal dist.']

for name, func, exact in zip(list_names,list_functions,list_exact):

    list_integral = complex(0.,0.)*N.zeros_like(list_hw)

    for i,wedge in enumerate(grid.list_wedges):

        integrand = func(wedge.list_k)

        list_Fk =  integrand[:,N.newaxis]*list_hw[N.newaxis,:]

        list_integral += AreaIntegrator(wedge,list_Fk)


    integral = N.real(N.average(list_integral/list_hw))

    error = N.abs(exact-integral)


    print ' %12s    %10.6e      %10.6e   %10.6e'%(name,exact,integral,error)


