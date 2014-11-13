#================================================================================
#
#           Module Compute_Loop_Function
#           ==========================================
#
#       This module computes the loop function, H^{alpha,nu}(q,hw), for a given nu,q.
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

    def __init__(self, mu, beta, q_vector, E_phonon_polarization, hw_ph, grid, external_list_hw, Gamma_width):

        self.mu      = mu
        self.beta    = beta

        self.q       = deepcopy(q_vector)
        self.E_ph    = deepcopy(E_phonon_polarization)
        self.hw_ph   = deepcopy(hw_ph)
        self.grid    = deepcopy(grid)

        self.list_hw = external_list_hw 
        self.nhw     = len(external_list_hw)

        # Widths
        self.Gamma   = Gamma_width

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
        # second index  : elph coupling parameters u, v
        self.Hq = complex(0.,0.)*N.zeros([2, 2, self.nhw])
        return


    def Compute_Hq(self):

        # Loop on wedges in the 1BZ
        for wedge in self.grid.list_wedges:

            # Compute matrix elements
            Mq = MGridFunction(q_vector=self.q,E_phonon_polarization=self.E_ph,hw_nu_q=self.hw_ph,wedge=wedge)

            # Compute the I function, which contains frequency dependence
            Iq = IGridFunction(q_vector=self.q,list_hw=self.list_hw, delta_width=self.Gamma, mu=self.mu, beta=self.beta, wedge=wedge)

            for i_alpha, alpha in zip([0,1],['x','y']):
                for i in [0,1]:
                    for n1 in [0,1]:
                        for n2 in [0,1]:
                            for n3 in [0,1]:

                                M_key = (alpha,i,n1,n2,n3)
                                M_index = Mq.index_dictionary[M_key]

                                I_key = (i,n1,n2,n3)
                                I_index = Iq.index_dictionary[I_key]

                                # dimension nk
                                MatrixElements_u = Mq.M[M_index,0,:]
                                MatrixElements_v = Mq.M[M_index,1,:]

                                # dimensions [nk,nw]
                                IElements = Iq.I[I_index,:,:]

                                MI_u = MatrixElements_u[:,N.newaxis]*IElements 
                                MI_v = MatrixElements_v[:,N.newaxis]*IElements 

                                self.Hq[i_alpha,0,:] += self.normalization*AreaIntegrator(wedge,MI_u)
                                self.Hq[i_alpha,1,:] += self.normalization*AreaIntegrator(wedge,MI_v)

        return
