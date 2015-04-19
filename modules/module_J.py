#================================================================================
#
#           Module J
#           ===============
#
#       This module implements an object which represents the frequency 
#       dependent term J_{ {n}} (k, k+q, w).
#
#================================================================================

from module_Constants import *
from module_Grids import *
from module_Functions import *
from module_HardCoreKernel  import *

import scipy.special as SS


class JGridFunction:
    """
    Class which builds and contains the J_{{n}}(k,k+q,w) object, which
    is the result of summing on all the G functions in the loop function.

    All the band energies will be assumed to be in eV, and the gamma 
    functions are in eV x a0^2. The units of J are
    
        [ J ] ~ Length^2/ Energy

    As computed, we'll have [ J ] ~ a0^2/eV. To convert to a0^2/Ha (atomic units),
    we must multipy by conversion factor = Ha_to_eV
    """

    def __init__(self, q_vector, hw, Green_Gamma_width, kernel_Gamma_width, mu, beta, wedge):

        self.conversion_factor = Ha_to_eV

        self.q_vector    = q_vector
        self.hw          = hw


        self.Green_Gamma_width  = Green_Gamma_width
        self.kernel_Gamma_width = kernel_Gamma_width

        # Whenever the imaginary part in the denominator vanishes,
        #  we'll use a Fadeeva broadening instead of a lorentzian broadening.
        self.fadeeva_delta_width = 0.025 # eV, to avoid divergences

        self.mu  = mu
        self.beta= beta

        # Create the scattering kernel object
        self.SK = ScatteringKernel(self.mu,self.beta,self.kernel_Gamma_width,self.Green_Gamma_width)

        # let's assume a default filename; read in the needed parameters
        filename ='scattering_spline.nc'
        self.SK.read_spline_parameters(filename)

        self.build_indices()

        self.nk   = len(wedge.list_k)
        self.dim  = len(self.index_dictionary)

        # J will be in units of [Length]^2 / [Enegy].
        # where all terms are in fundamental units

        self.J    = complex(0.,0.)*N.zeros([self.dim,self.nk])

        self.build_J(wedge)

    def build_indices(self):

        self.index_dictionary = {}
                
        index = -1
        for n1 in [0,1]:
            for n2 in [0,1]:
                for n3 in [0,1]:
                    index += 1 
                    key = (n1,n2,n3)
                    self.index_dictionary[key] = index

    def get_D(self,list_xi, list_xi2):
        """
        This routine computes the "D" function, which is fairly complicated to reproduce here.

        We have

        J_{n}(k,k+q,w) =                 D(xi_{n1 k},eta hbar omega)-D(xi_{n3 k+q}, hbar omega)
                                         ---------------------------------------------------------
                                                        xi_{n1 k} - xi_{n3 k+q}

        We'll assume that list_xi  and list_xi2 have the same dimension [nk].

        The equations implemented here are obtained in my document computing_J.pdf .
        """

        # Initialize D
        D = complex(0.,0.)*N.zeros_like(list_xi)

        # Build some ingredients which we'll use over and over.
        xi1  = 1.*list_xi
        xi2  = 1.*list_xi2


        for eta in [-1,1]:

            KK_fKR_xi1  = self.SK.get_integral_KK('fKR', xi1, sign_Gamma= eta)
            KK_fKI_xi1  = self.SK.get_integral_KK('fKI', xi1, sign_Gamma= eta)

            KK_dfKR_xi1  = self.SK.get_integral_KK('dfKR', xi1, sign_Gamma= eta)
            KK_dfKI_xi1  = self.SK.get_integral_KK('dfKI', xi1, sign_Gamma= eta)


            for phi in [-1,1]:

                #============
                # First term
                #============

                real_denominator    = xi1-xi2+phi*self.hw
                one_on_denominator  = N.real(self.cutoff_denominator(real_denominator,fadeeva_width=1e-6))

                prefactor = 0.5*complex(1.,0.)*phi*one_on_denominator  


                D += prefactor* ( eta*KK_fKR_xi1   \
                                 -eta*self.SK.get_integral_KK('fKR', xi2-phi*self.hw, sign_Gamma= eta) \
                                 +1j*KK_fKI_xi1 \
                                 -1j*self.SK.get_integral_KK('fKI', xi2-phi*self.hw, sign_Gamma= eta) )

                #============
                # second term
                #============

                real_denominator    = xi1-xi2+eta*self.hw
                if eta == phi:
                    one_on_denominator  = N.real(self.cutoff_denominator(real_denominator,fadeeva_width=1e-6))
                else:
                    one_on_denominator = 1./(real_denominator + 1j*(eta-phi)*self.Green_Gamma_width)

                prefactor = 0.5*complex(1.,0.)*phi*self.hw*one_on_denominator  

                D += prefactor* ( KK_dfKR_xi1   \
                                 -self.SK.get_integral_KK('dfKR', xi2-eta*self.hw, sign_Gamma= phi) \
                                 +1j*eta*KK_dfKI_xi1 \
                                 -1j*eta*self.SK.get_integral_KK('dfKI', xi2-eta*self.hw, sign_Gamma= phi) )


        return D


    def cutoff_denominator(self,list_x,fadeeva_width):
        """
        This function returns an approximation to 1/(x+i delta) using the Faddeeva function
        """
        z   = list_x/fadeeva_width

        den = -1j*SS.wofz(z)/fadeeva_width*N.sqrt(N.pi)

        return den


    def build_J(self,wedge):
        """
        All quantities are built here. It will be *crucial* to make
        this code as transparent as possible, to avoid making mistakes!
        """

        list_k    = wedge.list_k
        list_kq   = list_k + self.q_vector

        eps_k     = function_epsilon_k(list_k)

        list_epsilon_k = [-eps_k, eps_k]

        eps_kq    = function_epsilon_k(list_kq)
        list_epsilon_kq = [-eps_kq, eps_kq]


        for n2, list_epsilon2 in zip([0,1],list_epsilon_k):

            list_xi2 = list_epsilon2 - self.mu 

            list_D1 = []
            list_D3 = []

            # Compute the C functions
            for n1, list_epsilon1 in zip([0,1],list_epsilon_k):
                list_xi1 = list_epsilon1 - self.mu
                D1       = self.get_D(list_xi1, list_xi2)

                list_D1.append(D1)

            for n3, list_epsilon3 in zip([0,1],list_epsilon_kq):
                list_xi3 = list_epsilon3 - self.mu
                D3       = self.get_D(list_xi3, list_xi2)
                list_D3.append(D3)


            for n1, list_epsilon1, D1 in zip([0,1],list_epsilon_k, list_D1):

                list_xi1 = list_epsilon1 - self.mu 

                for n3, list_epsilon3, D3 in zip([0,1],list_epsilon_kq,list_D3):

                    list_xi3 =  list_epsilon3 - self.mu

                    key   = (n1,n2,n3)
                    index = self.index_dictionary[key]

                    # cutoff the 1/0 
                    den13 = N.real(self.cutoff_denominator(list_xi1-list_xi3,fadeeva_width=1e-6))

                    contribution = den13*(D1-D3) 

                    # add the "smooth" part to the J function
                    self.J[index,:] += self.conversion_factor*contribution

        return 



