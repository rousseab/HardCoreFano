#================================================================================
#
#           Module Compute_Loop_Function_Product_weak_scattering
#           ======================================================
#
#       This module computes the product of loop functions, H^{alpha,nu}(q,hw), 
#       for a given nu,q, and averages over all indices. 
#
#       The computation assumes the weak scattering limit, where gamma(iwm) --> I0,
#       for I0 "weak".
#       
#================================================================================

from module_Constants import *
from module_Grids import *
from module_Functions import *
from module_M import *
from module_I import *
from module_HardCoreKernel import *
from module_Integrators import *

class Compute_Loop_Function_Product_weak_scattering:

    def __init__(self, type_of_integral, mu, beta, q_vector, E_phonon_polarization, hw_ph, grid, \
                    external_list_hw,  delta_width):

        self.type_of_integral = type_of_integral

        self.mu      = mu
        self.beta    = beta

        self.q       = deepcopy(q_vector)
        self.E_ph    = deepcopy(E_phonon_polarization)
        self.hw_ph   = deepcopy(hw_ph)

        self.grid    = deepcopy(grid)

        self.list_hw = external_list_hw 
        self.nhw     = len(external_list_hw)

        # Widths
        self.delta_width = delta_width


        # Normalization contains change of measure in going from sum to integral
        # on k, also sum on spin:
        #
        #   1/Omega \sum_{k sigma} (...) -->  2/Omega x Omega/(2 pi)^2 int d2k (...)
        #                                    --> 1/2 pi^2 int d2k (...)

        self.normalization = 1./(2.*N.pi**2)

        # columns loop faster in C (and thus in python)
        # First index is for x or y, second index is for frequency

        self.initialize_Hq()

    def initialize_Hq(self):
        """
        Hq will be in fundamental units of [charge] x [velocity], e hbar / m a_0
        """
        # Hq will be in fundamental units of [charge] x [velocity], e hbar / m a_0

        # first index   : alpha = x, y
        # second index  : L = 'A','B' the impurity scattering site index
        self.Hq_plus  = complex(0.,0.)*N.zeros([ 2, 2, self.nhw])
        self.Hq_minus = complex(0.,0.)*N.zeros([ 2, 2, self.nhw])

        return

    def Compute_Hq_Average(self):

        # compute the smooth grid contribution
        self.Compute_Hq_Product(self.grid)

        # Compute the averaged product of H(q,w) H(-q,w)
        self.HqHq = 0.25* ( self.Hq_plus[0,0, :]*self.Hq_minus[0,0, :] +
                            self.Hq_plus[1,0, :]*self.Hq_minus[1,0, :] +
                            self.Hq_plus[0,1, :]*self.Hq_minus[0,1, :] +
                            self.Hq_plus[1,1, :]*self.Hq_minus[1,1, :] )


    def Compute_Hq_Product(self,grid):

        # Loop on wedges in the 1BZ
        for wedge in grid.list_wedges:

            # Compute matrix elements
            Mq = MGridFunction(q_vector=self.q,E_phonon_polarization=self.E_ph,hw_nu_q=self.hw_ph,wedge=wedge)

            # Compute the I function, which contains frequency dependence
            Iq = IGridFunction_weak_scattering(q_vector=self.q, list_hw=self.list_hw, \
                                delta_width=self.delta_width, mu=self.mu, beta=self.beta,  wedge=wedge)

            for n1 in [0,1]:
                for n2 in [0,1]:
                    for n3 in [0,1]:

                        I_key = (n1,n2,n3)
                        I_index = Iq.index_dictionary[I_key]
                        # dimensions [nk,nw]
                        IElements = Iq.I[I_index,:,:]

                        for i_alpha, alpha in zip([0,1],['x','y']):
                            for i_L, L in zip([0,1],['A','B']):


                                M_key   = (alpha,L,n1,n2,n3)
                                M_index = Mq.index_dictionary[M_key]


                                # dimension nk
                                MatrixElements_u = Mq.M[M_index,0,:]
                                # MatrixElements_v = Mq.M[M_index,1,:] no longer in use

                                MatrixElements_u_star = N.conjugate(MatrixElements_u)

                                MI      = MatrixElements_u[:,N.newaxis]*IElements 
                                M_starI = MatrixElements_u_star[:,N.newaxis]*IElements 

                                self.Hq_plus[i_alpha,i_L, :]  += self.normalization*AreaIntegrator(wedge,MI)
                                self.Hq_minus[i_alpha,i_L, :] += -self.normalization*AreaIntegrator(wedge,M_starI)



        return


