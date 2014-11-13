#================================================================================
#
#           Module I
#           ===============
#
#       This module implements an object which represents the frequency 
#       dependent term I_{i, {n}} (k, k+q, w).
#
#================================================================================

from module_Constants import *
from module_Grids import *
from module_Functions import *
from module_HardCoreKernel  import *

import scipy.special as SS


class IGridFunction:
    """
    Class which builds and contains the I_{i,{n}}(k,k+q,w) object, which
        is the result of summing on all the G functions in the loop function.

    All the band energies will be assumed to be in eV, and the gamma 
    functions are in eV x a0^2. The units of I are
    
        [ I ] ~ Length^2/ Energy

    As computed, we'll have [ I ] ~ a0^2/eV. To convert to a0^2/Ha (atomic units),
    we must multipy by conversion factor = Ha_to_eV

    """

    def __init__(self,q_vector,list_hw, delta_width, mu, beta, wedge):


                self.conversion_factor = Ha_to_eV

                self.q_vector    = q_vector
                self.delta_width = delta_width
                self.list_hw     = list_hw
                self.z           = list_hw+1j*delta_width
                self.z2          = self.z**2

                self.mu  = mu
                self.beta= beta


                self.build_indices()

                self.nk   = len(wedge.list_k)
                self.nw   = len(self.list_hw)
                self.dim  = len(self.index_dictionary)

        # I will be in units of [Length]^2 / [Enegy].
        # where all terms are in fundamental units

                self.I    = complex(0.,0.)*N.zeros([self.dim,self.nk,self.nw])

                self.epsilon_k = complex(0.,0.)*N.zeros([2,self.nk])
                self.epsilon_kq= complex(0.,0.)*N.zeros([2,self.nk])

                self.build_I(wedge)


    def build_indices(self):

        self.index_dictionary = {}
                
        index = -1
        for i in [0,1]:
            for n1 in [0,1]:
                for n2 in [0,1]:
                    for n3 in [0,1]:
                        index += 1 
                        key = (i,n1,n2,n3)
                        self.index_dictionary[key] = index


    def get_frequency_term(self, list_energy_difference):

        de  = list_energy_difference[:,N.newaxis]
        de2 = de**2

        freq = -2.*self.z/(de2-self.z2)

        return freq

    def cutoff_denominator(self, list_energy):
        """
        Compute 1/x, with a gaussian cutoff as we approach x=0. 
        """

        tol = 1e-3

        z   = list_energy/tol

        den = N.imag(SS.wofz(z))/tol*N.sqrt(N.pi)

        return den


    def build_I(self,wedge):
        """
        All quantities are built here. It will be *crucial* to make
        this code as transparent as possible, to avoid making mistakes!
        """

        list_k    = wedge.list_k
        list_kq   = list_k + self.q_vector

        fk        = function_fk(list_k)
        eps_k     = function_epsilon_k(list_k)

        list_epsilon_k = [-eps_k, eps_k]

        fkq       = function_fk(list_kq)
        eps_kq    = function_epsilon_k(list_kq)

        list_epsilon_kq = [-eps_kq, eps_kq]

        list_Fk  = []
        list_Fkq = []
        for epsilon_k, epsilon_kq in zip(list_epsilon_k,list_epsilon_kq):

            list_Fk.append(function_fermi_occupation(epsilon_k,self.mu,self.beta))
            list_Fkq.append(function_fermi_occupation(epsilon_kq,self.mu,self.beta))


        for n1, e1, f1 in zip([0,1],list_epsilon_k, list_Fk):

            for n3, e3, f3 in zip([0,1],list_epsilon_kq, list_Fkq):
                den13 = self.cutoff_denominator(e1-e3)

                for n2, e2, f2 in zip([0,1],list_epsilon_k, list_Fk):

                    for i, get_gamma in zip([0,1], [get_gamma1, get_gamma2]):

                        g1 = get_gamma(e1)
                        g3 = get_gamma(e3)

                        key   = (i, n1,n2,n3)
                        index = self.index_dictionary[key]

                        freq12 = self.get_frequency_term(e1-e2)

                        freq32 = self.get_frequency_term(e3-e2)

                        fac12  = (f1-f2)*g1
                        fac32  = (f3-f2)*g3

                        numerator = fac12[:,N.newaxis]*freq12 -fac32[:,N.newaxis]*freq32

                        self.I[index,:,:] = self.conversion_factor*den13[:,N.newaxis]*numerator 

        return 
