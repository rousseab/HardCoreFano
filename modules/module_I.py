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

    def __init__(self, q_vector, hw, Green_Gamma_width, kernel_Gamma_width, mu, beta, wedge):

        self.conversion_factor = Ha_to_eV

        self.q_vector    = q_vector
        self.hw          = hw


        self.Green_Gamma_width  = Green_Gamma_width
        self.kernel_Gamma_width = kernel_Gamma_width
        self.delta_width        = 0.001 # eV, to avoid divergences

        self.mu  = mu
        self.beta= beta

        # Create the scattering kernel object
        self.SK = ScatteringKernel(self.mu,self.beta,self.kernel_Gamma_width)

        # let's assume a default filename; read in the needed parameters
        filename ='scattering_spline.nc'
        self.SK.read_spline_parameters(filename)

        self.build_indices()

        self.nk   = len(wedge.list_k)
        self.dim  = len(self.index_dictionary)

        # I will be in units of [Length]^2 / [Enegy].
        # where all terms are in fundamental units

        self.I    = complex(0.,0.)*N.zeros([self.dim,self.nk])

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

    def get_C(self,list_xi, list_xi2, eta):
        """
        This routine computes the "C" function, which is fairly complicated to reproduce here.

        We have

        I_{n}(k,k+q,w) = -sum_{eta} eta  C(xi_{n1 k},eta hbar omega)-C(xi_{n3 k+q},eta hbar omega)
                                         ---------------------------------------------------------
                                                        xi_{n1 k} - xi_{n3 k+q}

        We'll assume that list_xi1 and list_xi2 have the same dimension [nk].
        """

        # Initialize C
        C = complex(0.,0.)*N.zeros_like(list_xi)

        # Build some ingredients which we'll use over and over.

        xi1  = list_xi
        xi2  = list_xi2
        xi2_minus_eta_hw = list_xi2-eta*self.hw

        real_denominator = xi1-xi2_minus_eta_hw 

        Fermi_1          = function_fermi_occupation(xi1, 0., self.beta)
        Fermi_2          = function_fermi_occupation(xi2, 0., self.beta)
        Fermi_2_minus_ehw= function_fermi_occupation(xi2_minus_eta_hw, 0., self.beta)

        KR_xi1          = get_KR(xi1+self.mu,self.kernel_Gamma_width)
        KI_xi1          = get_KI(xi1+self.mu,self.kernel_Gamma_width)

        KR_xi2_minus_ehw= get_KR(xi2_minus_eta_hw+self.mu,self.kernel_Gamma_width)
        KI_xi2_minus_ehw= get_KI(xi2_minus_eta_hw+self.mu,self.kernel_Gamma_width)

        fKI_KK_xi1           = self.SK.get_fKI_KK(xi1)
        fKI_KK_xi2_minus_ehw = self.SK.get_fKI_KK(xi2_minus_eta_hw)

        # Add each term, one by one. This may not be efficient, but it must be transparent

        denominator = real_denominator  + 1j*eta*self.Green_Gamma_width 
        numerator   = Fermi_2*( KR_xi2_minus_ehw + 1j*eta*KI_xi2_minus_ehw)-Fermi_1*KR_xi1
        C +=  numerator/denominator


        denominator = real_denominator  + 1j*eta*self.delta_width  + 1j*(eta+1.)*self.Green_Gamma_width 
        numerator   = -0.5j*( Fermi_1*KI_xi1 + eta*Fermi_2_minus_ehw*KI_xi2_minus_ehw ) \
                      -0.5 *(fKI_KK_xi1 - fKI_KK_xi2_minus_ehw )
        C +=  numerator/denominator

        denominator = real_denominator  + 1j*eta*self.delta_width  + 1j*(eta-1.)*self.Green_Gamma_width 
        numerator   = -0.5j*(-Fermi_1*KI_xi1 + eta*Fermi_2_minus_ehw*KI_xi2_minus_ehw ) \
                      -0.5 *(fKI_KK_xi1 - fKI_KK_xi2_minus_ehw )
        C +=  numerator/denominator

        return C

    def cutoff_denominator(self,list_x,fadeeva_width):
        """
        This function returns an approximation to 1/(x+i delta) using the Faddeeva function
        """
        z   = list_x/fadeeva_width

        den = -1j*SS.wofz(z)/fadeeva_width*N.sqrt(N.pi)

        return den


    def build_I(self,wedge):
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

            for eta in [-1.,1.]:

                list_C1 = []
                list_C3 = []

                # Compute the C functions
                for n1, list_epsilon1 in zip([0,1],list_epsilon_k):
                    list_xi1 = list_epsilon1 - self.mu
                    C1       = self.get_C(list_xi1, list_xi2, eta)
                    list_C1.append(C1)

                for n3, list_epsilon3 in zip([0,1],list_epsilon_kq):
                    list_xi3 = list_epsilon3 - self.mu
                    C3       = self.get_C(list_xi3, list_xi2, eta)
                    list_C3.append(C3)


                for n1, list_epsilon1, C1 in zip([0,1],list_epsilon_k, list_C1):

                    list_xi1 = list_epsilon1 - self.mu 

                    for n3, list_epsilon3, C3 in zip([0,1],list_epsilon_kq,list_C3):

                        list_xi3 =  list_epsilon3 - self.mu


                        key   = (n1,n2,n3)
                        index = self.index_dictionary[key]

                        # cutoff the 1/0 
                        den13 = self.cutoff_denominator(list_xi1-list_xi3,fadeeva_width=1e-6)

                        contribution = -eta*den13*(C1-C3) 


                        # add the "smooth" part to the I function
                        self.I[index,:] += self.conversion_factor*contribution

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

