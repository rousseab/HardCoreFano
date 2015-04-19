#================================================================================
#
#           Module Compute_Loop_Function
#           ==========================================
#
#       This module computes the loop function projected on the phonon frequency, 
#               H^{alpha L}_{nu q}(hw_{nu q}), for a given nu,q.
#
#       This is then stored as the unitless parameters R^{alpha L}_{nu q} and I^{alpha L}_{nu q}
#
#================================================================================

from module_Constants import *
from module_Grids import *
from module_Functions import *
from module_M import *
from module_J import *
from module_HardCoreKernel import *
from module_Integrators import *

class Compute_Loop_Function:

    def __init__(self, mu, beta, q_vector, E_phonon_polarization, hw_ph, 
                    grid, kernel_Gamma_width, Green_Gamma_width):

        self.mu      = mu
        self.beta    = beta

        self.q       = deepcopy(q_vector)
        self.E_ph    = deepcopy(E_phonon_polarization)
        self.hw_ph   = deepcopy(hw_ph)

        self.grid    = deepcopy(grid)


        # Widths
        self.Green_Gamma_width  = Green_Gamma_width
        self.kernel_Gamma_width = kernel_Gamma_width 


        # Normalization contains change of measure in going from sum to integral
        # on k, also sum on spin:
        #
        #   1/Omega \sum_{k sigma} (...) -->  2/Omega x Omega/(2 pi)^2 int d2k (...)
        #                                    --> 1/2 pi^2 int d2k (...)

        self.normalization = 1./(2.*N.pi**2)

        # initialize parameters
        self.initialize_R_and_I()

    def initialize_R_and_I(self):
        """
        Hq will be in fundamental units of [charge] x [velocity], e hbar / m a_0
        """
        # Hq will be in fundamental units of [charge] x [velocity], e hbar / m a_0

        # first index   : alpha = x, y
        # second index  : L = 'A','B' the impurity scattering site index

        self.Rq          = complex(0.,0.)*N.zeros([ 2, 2])
        self.Iq          = complex(0.,0.)*N.zeros([ 2, 2])

        return

    def Compute_R_and_I(self):

        # Loop on wedges in the 1BZ
        for wedge in self.grid.list_wedges:

            # Compute matrix elements
            Mq = MGridFunction(q_vector=self.q,E_phonon_polarization=self.E_ph,hw_nu_q=self.hw_ph,wedge=wedge)

            # Compute the I function, which contains frequency dependence
            Iq = IGridFunction(q_vector = self.q, hw = self.hw_ph, \
                               Green_Gamma_width  = self.Green_Gamma_width, \
                               kernel_Gamma_width = self.kernel_Gamma_width, \
                               mu = self.mu, beta = self.beta,  wedge = wedge)

            for n1 in [0,1]:
                for n2 in [0,1]:
                    for n3 in [0,1]:

                        I_key = (n1,n2,n3)
                        I_index = Iq.index_dictionary[I_key]
                        # dimensions [nk]
                        IElements = Iq.I[I_index,:]

                        for i_alpha, alpha in zip([0,1],['x','y']):
                            for i_L, L in zip([0,1],['A','B']):


                                M_key   = (alpha,L,n1,n2,n3)
                                M_index = Mq.index_dictionary[M_key]


                                # dimension [nk]
                                MatrixElements_u = Mq.M[M_index,0,:]

                                MatrixElements_u_star = N.conjugate(MatrixElements_u)

                                MI      = MatrixElements_u*IElements 
                                M_starI = MatrixElements_u_star*IElements 

                                # We use a trick to extract R and I
                                #
                                # We have H(q,w)  = sum M (R+iI)
                                #         H(-q,w) = sum M^* (-R-iI)
                                # such that
                                #         H(q,w)+H^*(-q,w) = sum M I
                                #         -----------------             
                                #               2i
                                # and      
                                #         H(q,w)-H^*(-q,w) = sum M R
                                #         -----------------             
                                #               2 

                                # The AreaIntegrator routine expects an array with shape [nk,nw].
                                # To leverage existing code, we artificially give MI a new dimension.
                                Hq_plus  =  self.normalization*AreaIntegrator(wedge,MI[:,N.newaxis] )[0] # only take one frequency
                                Hq_minus = -self.normalization*AreaIntegrator(wedge,M_starI[:,N.newaxis])[0]

                                self.Rq[i_alpha,i_L] +=  0.5   *( Hq_plus -   N.conjugate(Hq_minus) )
                                self.Iq[i_alpha,i_L] += -0.5*1j*( Hq_plus +   N.conjugate(Hq_minus) )


        return
