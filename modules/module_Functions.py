#================================================================================
#
#           Module Functions
#           ================
#
#   This module implements various functions which are useful to the
#   graphene tight binding model.
#
#================================================================================
from module_Constants import *

def function_fk(list_k):
    """
    This is the function f(k), which is just the sum of exponentials phases. 
    """

    exponent1  = 1j*N.dot(list_k,tau1_minus_tau3)
    exponent2  = 1j*N.dot(list_k,tau2_minus_tau3)

    phase1     = N.exp(exponent1)
    phase2     = N.exp(exponent2)
    phase3     = complex(1.,0.)*N.ones(len(list_k))

    fk  = phase1+phase2+phase3

    return fk

def function_epsilon_k(list_k):
    """
    Yields the positive dispersion branch, in the same units as |gamma1|.
    """

    fk     = function_fk(list_k)

    abs_fk = N.sqrt(N.real(fk)**2+N.imag(fk)**2)

    eps    = complex(1.,0.)*N.abs(tight_binding_gamma1)*abs_fk

    return eps

def function_Vk(list_k):

    phase1 = N.exp(1j*N.dot(list_k,tau1_minus_tau3))
    phase2 = N.exp(1j*N.dot(list_k,tau2_minus_tau3))
    phase3 = complex(1.,0.)*N.ones_like(list_k[:,0])

    Vkx    = hat_tau1[0]*phase1 + hat_tau2[0]*phase2 + hat_tau3[0]*phase3
    Vky    = hat_tau1[1]*phase1 + hat_tau2[1]*phase2 + hat_tau3[1]*phase3

    return Vkx, Vky

def function_fermi_occupation(eps_k,mu,beta):

    z   = beta*(eps_k- mu)

    zmin =-28. # Exp[zmin] ~ 10^{-12}
    zmax = 28. # Exp[zmax] ~ 10^12

    I_big   = z >= zmax
    I_rest  = (z < zmax)*(z > zmin)

    f_occ         = N.ones_like(z)
    f_occ[I_big]  = 0.
    f_occ[I_rest] = 1./(1.+N.exp(z[I_rest]))

    return complex(1.,0.)*f_occ

def d_Fermi_dxi(eps_k,mu,beta):

    z   = beta*(eps_k- mu)

    zmin =-28. # Exp[zmin] ~ 10^{-12}
    zmax = 28. # Exp[zmax] ~ 10^12

    I_big   = z >= zmax
    I_rest  = (z < zmax)*(z > zmin)

    df_dxi  = N.zeros_like(z)
    e       = N.exp(z[I_rest])

    df_dxi[I_rest] = -beta*e/(1.+e)**2

    return complex(1.,0.)*df_dxi

