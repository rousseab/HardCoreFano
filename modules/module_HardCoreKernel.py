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


gamma1_constant =  -2.*N.pi*hvF**2/D_cutoff
gamma2_constant =  N.sqrt(3.)*acell*D_cutoff/(12.*hvF)


def get_gamma1(list_epsilon):
    """
    The argument should be epsilon_{nk} = xi_{nk}+mu, in eV.

    The returned value is in eV x a0^2
    """
    gamma1 = complex(1.,0.)*gamma1_constant*list_epsilon/D_cutoff

    return gamma1

def get_gamma2(list_epsilon):
    """
    The argument should be epsilon_{nk} = xi_{nk}+mu, in eV

    The returned value is in eV x a0^2
    """

    gamma2 = complex(1.,0.)*gamma1_constant*gamma2_constant*N.ones_like(list_epsilon)

    return gamma2

