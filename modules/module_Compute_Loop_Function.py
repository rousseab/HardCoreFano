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

    def __init__(self, type_of_integral, mu, beta, q_vector, E_phonon_polarization, hw_ph, grid, external_list_hw, Gamma_width):

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
        # second index  : L = 'A','B' the impurity scattering site index
        # third index   : elph coupling parameters u, v
        self.Hq = complex(0.,0.)*N.zeros([2, 2, 2, self.nhw])
        return

    def Compute_Hq(self):

        # Loop on wedges in the 1BZ
        for wedge in self.grid.list_wedges:

            # Compute matrix elements
            Mq = MGridFunction(q_vector=self.q,E_phonon_polarization=self.E_ph,hw_nu_q=self.hw_ph,wedge=wedge)

            # Compute the I function, which contains frequency dependence
            Iq = IGridFunction(type_of_integral=self.type_of_integral, q_vector=self.q,list_hw=self.list_hw, \
                                delta_width=self.Gamma, mu=self.mu, beta=self.beta,  wedge=wedge)

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
                                MatrixElements_v = Mq.M[M_index,1,:]

                                MI_u = MatrixElements_u[:,N.newaxis]*IElements 
                                MI_v = MatrixElements_v[:,N.newaxis]*IElements 

                                self.Hq[i_alpha,i_L, 0,:] += self.normalization*AreaIntegrator(wedge,MI_u)
                                self.Hq[i_alpha,i_L, 1,:] += self.normalization*AreaIntegrator(wedge,MI_v)

        return

class Compute_Imaginary_Loop_Function:
    """
    This object computes the imaginary part of the loop function, which involves an integral
    on a product of two delta functions. This function is very difficult to sample numerically;
    here we perform the integral analytically in the cone approximation.

    Hq will be in fundamental units of [charge] x [velocity], e hbar / m a_0
    """
    def __init__(self, mu, beta, q_vector, E_phonon_polarization, hw_ph, external_list_hw,Gamma_width):

        self.mu      = mu           # in eV
        self.beta    = beta         # in eV
        self.Gamma   = Gamma_width  # in eV

        self.q       = deepcopy(q_vector)
        self.E_ph    = deepcopy(E_phonon_polarization)
        self.hw_ph   = deepcopy(hw_ph)  # in eV

        self.list_hw = external_list_hw # in eV
        self.nhw     = len(external_list_hw)

        self.initialize_imaginary_Hq()

        self.prepare_special_points()

        tol = 1e-8
        if N.linalg.norm(self.q) < tol:
            self.q_is_zero = True
        else:
            self.q_is_zero = False

            self.identify_Q()

        # units 1/(a0^3 eV^2)/Ha
        self.normalization = 1j/(2.*N.pi*hvF**3)*Ha_to_eV


    def initialize_imaginary_Hq(self):
        """
        Hq will be in fundamental units of [charge] x [velocity], e hbar / m a_0
        """
        # Hq will be in fundamental units of [charge] x [velocity], e hbar / m a_0

        # first index   : alpha = x, y
        # second index  : L = 'A','B' the impurity scattering site index
        # third index   : elph coupling parameters u, v
        self.Hq = complex(0.,0.)*N.zeros([2, 2, 2, self.nhw])
        return

    def prepare_special_points(self):
        """
        This routine tabulates Gamma and all the images of K1 and K2 in the 1BZ.
        """
        c = N.cos(N.pi/3.)
        s = N.sin(N.pi/3.)

        self.K1    = twopia*2./3.*N.array([ 1., 0.])
        self.K2    = twopia*2./3.*N.array([-1., 0.])

        self.list_special_points = twopia*2./3.*N.array([[ 0., 0.],  # Gamma
                                                         [ 1., 0.],  # K1
                                                         [-c , s ],  # 
                                                         [-c ,-s ],  # 
                                                         [-1., 0.],  # K2
                                                         [ c ,-s ],  # 
                                                         [ c , s ]]) # 

        self.Gamma_set = [0]
        self.K1_set    = [1,2,3]
        self.K2_set    = [4,5,6]

        return       

    def identify_Q(self):
        """
        Assume that q = Gamma + Q, or q = K + Q. This routine identifies the closest special
        vector, be it Gamma, K1 or K2.
        """

        Kmq = self.list_special_points - self.q

        # find index of smallest cartesian distance
        I = N.argmin( N.sum(Kmq**2,axis=1) )

        self.Q  = Kmq[I]
        self.nQ = N.linalg.norm(self.Q)


        self.Qhat = self.Q/N.linalg.norm(self.Q)

        z = N.array([0.,0.,1.])
        self.zQhat = N.cross(z,self.Qhat)[:2]

        self.q_is_Gamma = False
        self.q_is_K1    = False
        self.q_is_K2    = False

        # identify which special point q belongs to.
        # Also identify around which cones k will be so that both k and k+q will be near the 
        # Fermi surface.
        if I in self.Gamma_set:
            self.q_is_Gamma = True
            self.list_K_for_k = N.array([self.K1,self.K2])
        elif I in self.K1_set:
            self.q_is_K1    = True
            self.list_K_for_k = N.array([self.K1])
        elif I in self.K2_set:
            self.q_is_K2    = True
            self.list_K_for_k = N.array([self.K2])

    def get_cos_theta(self,eta_hw):
        """
        Compute cos(theta(Q,w,eta)), with the understanding that when this 
        number is not between -1 and 1, then there is no solution and 
        the integral should vanish.

        This function also accounts for the fact that xi_{n3 k+q} must be of
        the same sign as mu-eta hw, which is imposed by delta( epsilon_{n3 k+q}-mu +eta hw).

        We'll assume that sign[mu -eta hw] < 0 implies lower band, and 
                          sign[mu -eta hw] > 0 implies upper band.
        """

        k_mu     = N.abs(self.mu)/hvF
        k_ew     = N.abs(self.mu-eta_hw)/hvF

        list_cos_theta = ( k_ew**2-k_mu**2-self.nQ**2 )/(2.*k_mu*self.nQ)

        list_sin_theta = N.zeros_like(list_cos_theta) 

        bound_cos_condition = (list_cos_theta  > -1.0)*(list_cos_theta  <  1.0)

        upper_band_condition= self.mu-eta_hw >= 0.
        lower_band_condition= self.mu-eta_hw <  0.

        I_solutions            = N.where (  bound_cos_condition )[0]
        I_solutions_upper_band = N.where (  bound_cos_condition*upper_band_condition )[0]
        I_solutions_lower_band = N.where (  bound_cos_condition*lower_band_condition )[0]

        list_sin_theta[I_solutions] = N.sqrt(1.-list_cos_theta[I_solutions]**2)

        return I_solutions,I_solutions_upper_band,I_solutions_lower_band,   list_cos_theta, list_sin_theta


    def cutoff_denominator(self,list_x,fadeeva_width):
        """
        This function returns an approximation to 1/(x+i delta) using the Faddeeva function
        """
        z   = list_x/fadeeva_width

        den = complex(1.,0.)*N.real(-1j*SS.wofz(z)/fadeeva_width*N.sqrt(N.pi))

        return den

    def get_frequency_prefactor(self,eta_hw,I_solutions,list_sin_theta):
        """
        The frequency prefactor has units of  a0^3 eV^2

        [ mu -eta hw] ~ eV
        [  KR  ]      ~ a0^2 eV
        [  Q  ]       ~ 1/a0

        """

        frequency_prefactor  = N.zeros_like(eta_hw)

        frequency_prefactor[I_solutions] = N.abs((self.mu-eta_hw[I_solutions])/list_sin_theta[I_solutions]/self.nQ)\
                                            *N.real( get_gamma(self.mu-eta_hw[I_solutions]+1j*self.Gamma) ) 


        return frequency_prefactor


    def build_list_K(self,I_solutions,list_cos_theta,list_sin_theta):


        k_mu = N.abs(self.mu)/hvF

        c = k_mu*list_cos_theta[I_solutions][:,N.newaxis]*self.Qhat
        s = k_mu*list_sin_theta[I_solutions][:,N.newaxis]*self.zQhat

        list_K_plus  = c+s
        list_K_minus = c-s


        return list_K_plus, list_K_minus

    def Compute_Hq(self):
        """
        Compute the imaginary part of Hq, within the cone approximation
        """

        # find the index n2, which is known immediately from the sign of the chemical potential 

        if self.q_is_zero:
            # the result is assumed to be zero 
            return

        if self.mu >= 0.:
            n2 = 1
        else:            
            n2 = 0

        for eta in [-1.,1.]:
            eta_hw = eta*self.list_hw                


            frequency_prefactor = complex(0.,0.)*N.zeros_like(eta_hw)


            I_solutions, I_solutions_upper_band, \
            I_solutions_lower_band, list_cos_theta, list_sin_theta = self.get_cos_theta(eta_hw)

            if len(I_solutions) == 0:
                # no contributions for this eta! Move on  
                continue

            # extract the kernel function and delta-function derived metric
            frequency_prefactor = self.get_frequency_prefactor(eta_hw,I_solutions,list_sin_theta)

            for Isol, n3 in zip([I_solutions_lower_band, I_solutions_upper_band],[0,1]):
                # Build the K vectors where the double deltas are peaked
                list_K_plus, list_K_minus =  self.build_list_K(Isol,list_cos_theta,list_sin_theta)

                for special_K in self.list_K_for_k:

                    list_k_plus = list_K_plus+special_K 
                    list_k_minus= list_K_minus+special_K 

                    DummyWedgePlus = DummyWedge(list_k_plus )
                    DummyWedgeMinus= DummyWedge(list_k_minus)

                    # build the needed matrix elements

                    MqPlus = MGridFunction(q_vector=self.q,E_phonon_polarization=self.E_ph,\
                                        hw_nu_q=self.hw_ph,wedge=DummyWedgePlus)

                    MqMinus= MGridFunction(q_vector=self.q,E_phonon_polarization=self.E_ph,\
                                        hw_nu_q=self.hw_ph,wedge=DummyWedgeMinus)


                    for n1 in [0,1]:

                        if n1 == n2:
                            xi1_factor = 0.
                        else:
                            xi1_factor = -2.*self.mu
                        
                        # Careful! do not divide by zero!!!

                        # units: [ frequency_prefactor] ~ a0^3 eV^2 
                        F = frequency_prefactor*eta_hw*self.cutoff_denominator(xi1_factor+eta_hw,fadeeva_width=0.001)

                        for i_alpha, alpha in zip([0,1],['x','y']):
                            for i_L, L in zip([0,1],['A','B']):

                                M_key   = (alpha,L,n1,n2,n3)
                                M_index = MqPlus.index_dictionary[M_key]

                                # dimension  Isol
                                MatrixElements_u = MqPlus.M[M_index,0,:]+MqMinus.M[M_index,0,:]
                                MatrixElements_v = MqPlus.M[M_index,1,:]+MqMinus.M[M_index,1,:]

                                # units: 
                                #        Matrix elements: ~ e hbar^3/(m^2 a0^3)
                                #        normalization  : ~ 1/(a0^3 eV^2)/Ha
                                #               F       : ~ a0^3 eV^2 
                                #               H       : ~ e hbar/ (m a0)
                                self.Hq[i_alpha,i_L, 0,Isol] += self.normalization*MatrixElements_u*F[Isol]
                                self.Hq[i_alpha,i_L, 1,Isol] += self.normalization*MatrixElements_v*F[Isol]

        return

class Compute_Loop_Function_MATSUBARA_SUMS:

    def __init__(self, mu, beta, q_vector, E_phonon_polarization, hw_ph, grid, external_list_hw, Gamma_width, matsubara_cutoff_energy):

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


        # cutoff for the generation of the matrsubara frequencies
        self.matsubara_cutoff_energy = matsubara_cutoff_energy

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
        # third index   : elph coupling parameters u, v
        self.Hq = complex(0.,0.)*N.zeros([2, 2, 2, self.nhw])
        return


    def Compute_Hq(self):

        # Loop on wedges in the 1BZ
        for wedge in self.grid.list_wedges:

            # Compute matrix elements
            Mq = MGridFunction(q_vector=self.q,E_phonon_polarization=self.E_ph,hw_nu_q=self.hw_ph,wedge=wedge)

            # Compute the I function, which contains frequency dependence
            Iq = ICGridFunction(q_vector=self.q,list_hw=self.list_hw, delta_width=self.Gamma, mu=self.mu,  \
                                    beta=self.beta, matsubara_grid_energy_cutoff=self.matsubara_cutoff_energy, wedge=wedge)

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
                                MatrixElements_v = Mq.M[M_index,1,:]

                                MI_u = MatrixElements_u[:,N.newaxis]*IElements 
                                MI_v = MatrixElements_v[:,N.newaxis]*IElements 

                                self.Hq[i_alpha,i_L, 0,:] += self.normalization*AreaIntegrator(wedge,MI_u)
                                self.Hq[i_alpha,i_L, 1,:] += self.normalization*AreaIntegrator(wedge,MI_v)

        return


