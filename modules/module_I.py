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
    Class which builds and contains the I_{{n}}(k,k+q,w) object, which
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
        for n1 in [0,1]:
            for n2 in [0,1]:
                for n3 in [0,1]:
                    index += 1 
                    key = (n1,n2,n3)
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

                    g1 = get_gamma(e1)
                    g3 = get_gamma(e3)

                    key   = (n1,n2,n3)
                    index = self.index_dictionary[key]

                    freq12 = self.get_frequency_term(e1-e2)

                    freq32 = self.get_frequency_term(e3-e2)

                    fac12  = (f1-f2)*g1
                    fac32  = (f3-f2)*g3

                    numerator = fac12[:,N.newaxis]*freq12 -fac32[:,N.newaxis]*freq32

                    self.I[index,:,:] = self.conversion_factor*den13[:,N.newaxis]*numerator 

        return 


class ICGridFunction:
    """
    Class which builds and contains the IC_{{n}}(k,k+q,w) object, which
    is the result of summing on all the G functions in the loop function,
    performing explicitly the Matsubara sum over fermionic frequencies.

    All the band energies will be assumed to be in eV, and the gammaC
    function is in eV x a0^2. The units of IC are
    
        [ IC ] ~ Length^2/ Energy

    As computed, we'll have [ IC ] ~ a0^2/eV. To convert to a0^2/Ha (atomic units),
    we must multipy by conversion factor = Ha_to_eV
    """

    def __init__(self,q_vector,list_hw, delta_width, mu, beta, matsubara_grid_energy_cutoff, wedge):


        self.conversion_factor = Ha_to_eV

        self.q_vector    = q_vector
        self.delta_width = delta_width
        self.list_hw     = list_hw
        self.z           = list_hw+1j*delta_width
        self.z2          = self.z**2

        self.mu  = mu
        self.beta= beta

        # matsubara frequencies  are of the form wm = pi/beta(2m+1)
        # let's find the largest m for the given cutoff
        m_max = N.int(N.ceil(beta*matsubara_grid_energy_cutoff/(2.*N.pi)))

        # Define the Matsubara frequencies, symmetric about 0, all the way to the cutoff
        self.iwm = 1j*N.pi/self.beta*(2*N.arange(-m_max,m_max)+1)


        # build, once and for all, the values of the kernel at the matsubara frequencies
        self.gammaC = get_gammaC(self.iwm+self.mu)


        self.build_indices()

        self.nk   = len(wedge.list_k)
        self.nw   = len(self.list_hw)
        self.dim  = len(self.index_dictionary)

        # I will be in units of [Length]^2 / [Enegy].
        # where all terms are in fundamental units

        # Dimensions:  [ keys,  nk,  hw ]  
        self.I    = complex(0.,0.)*N.zeros([self.dim,self.nk,self.nw])

        self.epsilon_k = complex(0.,0.)*N.zeros([2,self.nk])
        self.epsilon_kq= complex(0.,0.)*N.zeros([2,self.nk])

        self.build_I(wedge)


    def build_indices(self):

        self.index_dictionary = {}
                
        index = -1
        for n1 in [0,1]:
            for n2 in [0,1]:
                for n3 in [0,1]:
                    index += 1 
                    key = (n1,n2,n3)
                    self.index_dictionary[key] = index


    def get_green_function_denominator(self, list_xi):
        """
        Produce a matrix of the form

        G[nk, iwm] = 1/(iwm- xi_k)

        """
        G = complex(1.,0.)/(list_xi[:,N.newaxis]- self.gammaC[N.newaxis,:] )

        return G

    def get_frequency_double_green_function(self, list_xi):
        """
        Produce the tensor
        D [ nk, hw, iwm]
        D = 1./( (iwm-xik)**2 - (hw+idelta)^2 ) 
        """

        t1 = (self.iwm[N.newaxis,N.newaxis,:]-list_xi[:,N.newaxis,N.newaxis])**2
        t2 = self.z2[N.newaxis,:,N.newaxis]

        D  = 1./(t1-t2)

        return D


    def build_I(self,wedge):
        """
        All quantities are built here. It will be *crucial* to make
        this code as transparent as possible, to avoid making mistakes!

        Note that numpy arrays are by default in row-major order. That means
        that for A[i,j,k], k loops fastest.
        """

        list_k    = wedge.list_k
        list_kq   = list_k + self.q_vector

        eps_k     = function_epsilon_k(list_k)

        list_epsilon_k = [-eps_k, eps_k]

        eps_kq    = function_epsilon_k(list_kq)
        list_epsilon_kq = [-eps_kq, eps_kq]

        # Compute the green's function 

        list_gammaC_Gk = []
        list_Gkq       = []
        for ek in list_epsilon_k:
            # The green function is a function of xi=e-mu                

            # G[nk, iwm] = 1/(iwm- xi_k)
            Gk = self.get_green_function_denominator(ek-self.mu) 
            list_gammaC_Gk.append( Gk[:,:]*self.gammaC[N.newaxis,:] )

        for ekq in list_epsilon_kq:
            # The green function is a function of xi=e-mu                
            list_Gkq.append( self.get_green_function_denominator(ekq-self.mu) )

        
        # Combine all the frequency dependent terms
        for n2, ek2 in zip([0,1],list_epsilon_k):

            # Dimensions [ nk, hw, iwm]
            frequency_double_G = self.get_frequency_double_green_function(ek2-self.mu)

            for n1, gammaC_Gk in zip([0,1],list_gammaC_Gk):
                for n3, Gkq in zip([0,1],list_Gkq):

                    key   = (n1,n2,n3)
                    index = self.index_dictionary[key]


                    # dimension [nk, iwm] 
                    green_product = gammaC_Gk[:,:]*Gkq[:,:]

                    # dimension [ nk, hw, iwm]
                    summand = frequency_double_G[:,:,:]*green_product[:,N.newaxis,:]

                    # dimension [ nk, hw], units a0^2/Ha
                    self.I[index,:,:] = -2./self.beta*self.conversion_factor*N.sum(summand,axis=2)

        return 


def test_new_axis_slicing():
    """
    We do some pretty funky array slicing in the code generating I, and in particular
    we use N.newaxis quite liberally. It will be crucial to make sure that the code does
    exactly what we think it does. The following little test basically confirms 
    that.

    We could eventually turn this test into a unit test.
    """

    # fill these arrays with random  numbers. Make sure they have different
    # dimensions
    list_iwm = N.random.random(10)
    list_gam = N.random.random(10)
    list_xi  = N.random.random(20)
    list_z2  = N.random.random(30)

    t1 = (list_iwm[N.newaxis,N.newaxis,:]-list_xi[:,N.newaxis,N.newaxis])**2
    t2 = list_z2[N.newaxis,:,N.newaxis]

    # dimension [xi,z2,iwm]
    D  = 1./(t1-t2)


    G = 1./(list_xi[:,N.newaxis]- list_iwm[N.newaxis,:] )

    # dimension [xi,iwm]
    gG = G[:,:]*list_gam[N.newaxis,:] 

    sum = N.sum( D[:,:,:]*gG[:,N.newaxis,:],axis=2)

    error1 = 0.
    error2 = 0.
    error3 = 0.
    error4 = 0.

    for i,xi in enumerate(list_xi):
        for j,z2 in enumerate(list_z2):
            for k,iwm in enumerate(list_iwm):

                explicit_D = 1./( (iwm-xi)**2 - z2)

                difference = explicit_D-D[i,j,k] 

                error1+= difference**2 

    for i,xi in enumerate(list_xi):
        for j,iwm in enumerate(list_iwm):
            explicit_G = 1./(xi-iwm)

            difference = explicit_G-G[i,j] 
            error2+= difference**2 

    for i,xi in enumerate(list_xi):
        for j,(iwm,gam) in enumerate(zip(list_iwm,list_gam)):

            explicit_G  = 1./(xi-iwm)
            explicit_gG = explicit_G*gam  

            difference = explicit_gG-gG[i,j] 
            error3+= difference**2 

    for i,xi in enumerate(list_xi):
        for j,z2 in enumerate(list_z2):

            explicit_sum = 0.
            for k,(iwm,gam) in enumerate(zip(list_iwm,list_gam)):

                explicit_D  = 1./( (iwm-xi)**2 - z2)
                explicit_G  = 1./(xi-iwm)
                explicit_gG = explicit_G*gam  

                explicit_sum += explicit_D*explicit_gG 


            difference = explicit_sum-sum[i,j] 

            error4 += difference**2 



    error1= N.sqrt(error1)
    error2= N.sqrt(error2)
    error3= N.sqrt(error3)
    error4= N.sqrt(error4)

    print 'ERROR 1: %12.4e'%error1
    print 'ERROR 2: %12.4e'%error2
    print 'ERROR 3: %12.4e'%error3
    print 'ERROR 4: %12.4e'%error4

    return

if __name__=='__main__':
    test_new_axis_slicing()

