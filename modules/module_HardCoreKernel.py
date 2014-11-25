#================================================================================
#
#           Module HardCoreKernel
#           =====================
#
#   This module implements the simplest analytical functions which describe
#   the hard core kernel for graphene.
#
#================================================================================
from module_Constants import *
from module_Functions import *


# Define a few useful constants which are only useful for the Hard Core Kernel

k_cutoff = N.sqrt(2.*N.pi/Area) # in 1/bohr
D_cutoff = hvF*k_cutoff         # in eV

gamma_constant =  -2.*N.pi*hvF**2/D_cutoff

def get_gamma(list_epsilon):
    """
    The argument should be epsilon_{nk} = xi_{nk}+mu, in eV.

    The returned value is in eV x a0^2
    """
    gamma = complex(1.,0.)*gamma_constant*list_epsilon/D_cutoff

    return gamma


def get_gammaA(list_z):
    """
    The "real" kernel function is given by

        gamma(z) =  2 pi (hvF)^2/D ( D/z 1/ln[1-D^2/z^2])

    This function is wildly non-analytic near zero. We thus wish to
    extract its behavior at inifnity.    

    z-> infinity limit:
    gamma_inf (z) =- 2 pi (hvF)^2/D  z/D

    Pade approximation of the leftover, to first order:
    gamma_1   (z) =  2 pi (hvF)^2/D  ( 3x^3+2x)/(6x^4+3x^2+2),  x  = z/D

    gammaA = gamma_inf + gamma_1 

    The argument of the function should always be offset by mu!

    The returned value is in eV x a0^2
    """

    prefactor =  complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff

    x = list_z/D_cutoff
    x2= x*x
    x3= x2*x
    x4= x3*x

    gamma_inf = -prefactor*x 

    numerator   = complex(3.,0.)*x3+complex(2.,0.)*x
    denominator = complex(6.,0.)*x4+complex(3.,0.)*x2+complex(2.,0.)

    gamma_1   = prefactor*numerator/denominator 

    gammaA = gamma_inf+ gamma_1   

    return gammaA

def get_gammaC(list_z):
    """
    The "real" kernel function is given by

        gamma(z) =  2 pi (hvF)^2/D ( D/z 1/ln[1-D^2/z^2])

    This function is wildly non-analytic near zero. 

    gammaC = gamma - gammaA

    The argument of the function should always be offset by mu!

    The returned value is in eV x a0^2
    """

    prefactor =  complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff

    gammaA = get_gammaA(list_z)

    x = complex(1.,0.)*list_z/D_cutoff
    x2= x*x

    denominator = x*N.log(complex(1.,0.)-complex(1.,0.)/x2)

    gamma = prefactor/denominator 

    gammaC = gamma - gammaA   

    return gammaC

