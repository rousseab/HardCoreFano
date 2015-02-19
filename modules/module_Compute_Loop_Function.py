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
from module_I import *
from module_HardCoreKernel import *
from module_Integrators import *

class Compute_Loop_Function:

    def __init__(self, mu, beta, q_vector, E_phonon_polarization, hw_ph, 
                    grid_singular, grid_smooth, kernel_Gamma_width, delta_width):

        self.mu      = mu
        self.beta    = beta

        self.q       = deepcopy(q_vector)
        self.E_ph    = deepcopy(E_phonon_polarization)
        self.hw_ph   = deepcopy(hw_ph)

        self.grid_singular    = deepcopy(grid_singular)
        self.grid_smooth      = deepcopy(grid_smooth)


        # we do this to leverage code which is already there!
        self.list_hw = N.array([hw_ph])


        # Widths
        self.delta_width = delta_width
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

        self.Rq_smooth   = complex(0.,0.)*N.zeros([ 2, 2])
        self.Iq_smooth   = complex(0.,0.)*N.zeros([ 2, 2])

        self.Rq_singular = complex(0.,0.)*N.zeros([ 2, 2])
        self.Iq_singular = complex(0.,0.)*N.zeros([ 2, 2])

        self.Rq          = complex(0.,0.)*N.zeros([ 2, 2])
        self.Iq          = complex(0.,0.)*N.zeros([ 2, 2])

        return

    def Compute_R_and_I(self):

        # compute the smooth grid contribution
        self.Compute_R_and_I_per_grid(self.grid_smooth,'smooth')

        # compute the singular grid contribution
        self.Compute_R_and_I_per_grid(self.grid_singular,'singular')

        # Compute the sum
        self.Rq = self.Rq_smooth + self.Rq_singular
        self.Iq = self.Iq_smooth + self.Iq_singular


    def Compute_R_and_I_per_grid(self,grid,type_of_integral):


        # Loop on wedges in the 1BZ
        for wedge in grid.list_wedges:

            # Compute matrix elements
            Mq = MGridFunction(q_vector=self.q,E_phonon_polarization=self.E_ph,hw_nu_q=self.hw_ph,wedge=wedge)

            # Compute the I function, which contains frequency dependence
            Iq = IGridFunction(type_of_integral=type_of_integral, q_vector=self.q,list_hw=self.list_hw, \
                                delta_width=self.delta_width,kernel_Gamma_width=self.kernel_Gamma_width, \
                                mu=self.mu, beta=self.beta,  wedge=wedge)

            for n1 in [0,1]:
                for n2 in [0,1]:
                    for n3 in [0,1]:

                        I_key = (n1,n2,n3)
                        I_index = Iq.index_dictionary[I_key]
                        # dimensions [nk,nw = 1]
                        IElements = Iq.I[I_index,:,0]

                        for i_alpha, alpha in zip([0,1],['x','y']):
                            for i_L, L in zip([0,1],['A','B']):


                                M_key   = (alpha,L,n1,n2,n3)
                                M_index = Mq.index_dictionary[M_key]


                                # dimension nk
                                MatrixElements_u = Mq.M[M_index,0,:]

                                MatrixElements_u_star = N.conjugate(MatrixElements_u)

                                MI      = MatrixElements_u[:,N.newaxis]*IElements 
                                M_starI = MatrixElements_u_star[:,N.newaxis]*IElements 

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


                                Hq_plus  =  self.normalization*AreaIntegrator(wedge,MI)[0] # only take one frequency
                                Hq_minus = -self.normalization*AreaIntegrator(wedge,M_starI)[0]

                                if type_of_integral == 'smooth':
                                    self.Rq_smooth[i_alpha,i_L] +=  0.5   *( Hq_plus -   N.conjugate(Hq_minus) )
                                    self.Iq_smooth[i_alpha,i_L] += -0.5*1j*( Hq_plus +   N.conjugate(Hq_minus) )

                                elif type_of_integral == 'singular':
                                    self.Rq_singular[i_alpha,i_L] +=  0.5   *( Hq_plus -   N.conjugate(Hq_minus) )
                                    self.Iq_singular[i_alpha,i_L] += -0.5*1j*( Hq_plus +   N.conjugate(Hq_minus) )

        return


