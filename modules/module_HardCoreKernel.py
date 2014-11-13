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

