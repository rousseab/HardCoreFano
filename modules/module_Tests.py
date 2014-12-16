#================================================================================
#
#           Module Tests
#           ===============
#
#       This module implements various testing objects.
#
#================================================================================

from module_Constants import *
from module_Grids import *
from module_Functions import *
from module_HardCoreKernel  import *

import scipy.special as SS


class dF_GridFunction:
    """
    Class which builds and contains the function

        dF = - ( f(xi_{k+q}) - f(xi_{k}) / (xi_{k+q} - xi_{k})

    which should look like a delta function for small q.

    """

    def __init__(self,q_vector, mu, beta, wedge):

        self.q_vector    = q_vector

        self.mu  = mu
        self.beta= beta

        self.build_dF(wedge)


    def cutoff_denominator(self, list_energy):
        """
        Compute 1/x, with a gaussian cutoff as we approach x=0. 
        """

        tol = 1e-3

        z   = list_energy/tol

        den = N.imag(SS.wofz(z))/tol*N.sqrt(N.pi)

        return den


    def build_dF(self,wedge):
        """
        The function -(f_{k+q}-f_{k})/(xi_k+q - xi_k)
        is computed here
        """

        list_k    = wedge.list_k
        list_kq   = list_k + self.q_vector

        fk        = function_fk(list_k)
        eps_k     = function_epsilon_k(list_k)

        Fk = function_fermi_occupation(-eps_k,self.mu,self.beta)

        fkq       = function_fk(list_kq)
        eps_kq    = function_epsilon_k(list_kq)

        Fkq = function_fermi_occupation(-eps_kq,self.mu,self.beta)

        den = self.cutoff_denominator(-eps_kq+eps_k)

        self.df   = -(Fkq-Fk)*den

        return 
