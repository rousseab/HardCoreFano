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

    def __init__(self, type_of_integral, q_vector,list_hw, delta_width, kernel_Gamma_width, mu, beta, wedge):


        if type_of_integral == 'smooth':
            self.smooth   = True
            self.singular = False
        elif type_of_integral == 'singular':
            self.smooth   = False
            self.singular = True
        else :
            print 'ERROR! Pick a type of integral: smooth or singular'

        self.conversion_factor = Ha_to_eV

        self.q_vector    = q_vector
        self.delta_width = delta_width
        self.list_hw     = list_hw


        self.kernel_Gamma_width = kernel_Gamma_width

        self.mu  = mu
        self.beta= beta

        # Create the scattering kernel object
        self.SK = ScatteringKernel(self.mu,self.beta,self.kernel_Gamma_width)

        # let's assume a default filename; read in the needed parameters
        filename ='scattering_spline.nc'
        self.SK.read_spline_parameters(filename)

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

    def cutoff_denominator(self,list_x,fadeeva_width):
        """
        This function returns an approximation to 1/(x+i delta) using the Faddeeva function
        """
        z   = list_x/fadeeva_width

        den = -1j*SS.wofz(z)/fadeeva_width*N.sqrt(N.pi)

        return den

    def get_Y_smooth(self,list_xi, list_xi2, list_eta_hw):
        """
        This routine computes the "smooth" part of the Y function,

        Y(xi,xi2,eta_hw) =  ( L(xi,0) -L(xi2-eta_hw, 0 ) ) / (xi-[xi2-eta_hw]).

        It is understood that as xi1 -> xi2-eta_hw, the function should remain smooth. 
        Care must be taken to avoid singularities.

        We'll assume that list_xi1 and list_xi2 have the same dimension [nk], and list_eta_hw is [nw]
        """

        L1  = self.SK.get_L_oneD(list_xi)[:,N.newaxis]

        u   = list_xi2[:,N.newaxis]-list_eta_hw[N.newaxis,:]

        shape = u.shape

        L2  =  self.SK.get_L_oneD(u.flatten()).reshape(shape)

        numerator = L1-L2

        denominator = list_xi[:,N.newaxis]-list_xi2[:,N.newaxis]+list_eta_hw[N.newaxis,:]


        # Where the denominator is too close to zero, replace by value for argument close to zero.
        tol = 1e-6

        I,J = N.nonzero( N.abs(denominator) < tol)
        
        # substitute with safe value, to not divide by zero        
        denominator[I,J] = 1.
        Y_smooth = numerator/denominator

        # replace by derivative 

        L1_safe_plus  = self.SK.get_L_oneD(list_xi[I]+tol)
        L1_safe_minus = self.SK.get_L_oneD(list_xi[I]-tol)

        # note: u[I,J] is a 1D array, not a 2D array. This is because I and J are arrays, not slicing operators
        
        L2_safe_plus  = self.SK.get_L_oneD(u[I,J]+tol)
        L2_safe_minus = self.SK.get_L_oneD(u[I,J]-tol)

        # I'm not sure why this strange derivative works better, but I find that it does for tol large (for debugging)
        Y_smooth[I,J] = 0.5*((L1_safe_plus+L2_safe_plus)-(L1_safe_minus+L2_safe_minus))/(2.*tol)

        return Y_smooth

    def build_I(self,wedge):
        if self.smooth:
            self.build_I_smooth(wedge)
            #self.build_I_smooth_intraband_only(wedge)
        elif self.singular:
            self.build_I_singular(wedge)
            #self.build_I_singular_intraband_only(wedge)

    def build_I_smooth(self,wedge):
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

            list_xi2 = N.real(  list_epsilon2 - self.mu )

            for eta in [-1.,1.]:

                list_eta_hw =  eta*self.list_hw

                list_Y12 = []
                list_Y32 = []

                for n1, list_epsilon1 in zip([0,1],list_epsilon_k):

                    list_xi1 = N.real(  list_epsilon1 - self.mu )

                    Y12 = self.get_Y_smooth(list_xi1, list_xi2, list_eta_hw)
                    list_Y12.append(Y12)

                for n3, list_epsilon3 in zip([0,1],list_epsilon_kq):

                    list_xi3 = N.real(  list_epsilon3 - self.mu )

                    Y32 = self.get_Y_smooth(list_xi3, list_xi2, list_eta_hw)
                    list_Y32.append(Y32)


                for n1, list_epsilon1, Y12 in zip([0,1],list_epsilon_k, list_Y12):

                    list_xi1 = N.real(  list_epsilon1 - self.mu )


                    for n3, list_epsilon3, Y32 in zip([0,1],list_epsilon_kq,list_Y32):

                        list_xi3 = N.real(  list_epsilon3 - self.mu )


                        key   = (n1,n2,n3)
                        index = self.index_dictionary[key]

                        # cutoff the 1/0 
                        den13 = N.real ( self.cutoff_denominator(list_xi1-list_xi3,fadeeva_width=1e-6) )

                        smooth_contribution = eta*den13[:,N.newaxis]*(Y12-Y32) 


                        # add the "smooth" part to the I function
                        self.I[index,:,:] += self.conversion_factor*smooth_contribution

        return 

    def build_I_smooth_intraband_only(self,wedge):
        """
        All quantities are built here. It will be *crucial* to make
        this code as transparent as possible, to avoid making mistakes!

        Only intraband contributions are considered, assuming mu < 0.
        """

        list_k    = wedge.list_k
        list_kq   = list_k + self.q_vector

        eps_k     = function_epsilon_k(list_k)

        list_epsilon_k = [-eps_k, eps_k]

        eps_kq    = function_epsilon_k(list_kq)
        list_epsilon_kq = [-eps_kq, eps_kq]


        for n2, list_epsilon2 in zip([0,1],list_epsilon_k):

            list_xi2 = N.real(  list_epsilon2 - self.mu )

            for eta in [-1.,1.]:

                list_eta_hw =  eta*self.list_hw

                list_Y12 = []
                list_Y32 = []

                for n1, list_epsilon1 in zip([0,1],list_epsilon_k):

                    list_xi1 = N.real(  list_epsilon1 - self.mu )

                    Y12 = self.get_Y_smooth(list_xi1, list_xi2, list_eta_hw)
                    list_Y12.append(Y12)

                for n3, list_epsilon3 in zip([0,1],list_epsilon_kq):

                    list_xi3 = N.real(  list_epsilon3 - self.mu )

                    Y32 = self.get_Y_smooth(list_xi3, list_xi2, list_eta_hw)
                    list_Y32.append(Y32)


                for n1, list_epsilon1, Y12 in zip([0,1],list_epsilon_k, list_Y12):

                    list_xi1 = N.real(  list_epsilon1 - self.mu )


                    for n3, list_epsilon3, Y32 in zip([0,1],list_epsilon_kq,list_Y32):

                        list_xi3 = N.real(  list_epsilon3 - self.mu )


                        key   = (n1,n2,n3)
                        index = self.index_dictionary[key]

                        # cutoff the 1/0 
                        den13 = N.real ( self.cutoff_denominator(list_xi1-list_xi3,fadeeva_width=1e-6) )

                        smooth_contribution = eta*den13[:,N.newaxis]*(Y12-Y32) 

                        if n1 == 0 and n2 == 0 and n3 == 0:
                            smooth_contribution_intraband = smooth_contribution 
                        else:
                            smooth_contribution_intraband = complex(0.,0.)*N.zeros_like(smooth_contribution)

                        # add the "smooth" part to the I function
                        self.I[index,:,:] += self.conversion_factor*smooth_contribution_intraband

        return 

    def build_I_singular(self,wedge):
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

            list_xi2 = N.real(  list_epsilon2 - self.mu )

            for eta in [-1.,1.]:

                list_eta_hw =  eta*self.list_hw

                KR = get_KR(list_xi2[:,N.newaxis]+self.mu-list_eta_hw[N.newaxis,:],self.kernel_Gamma_width)
                KI = get_KI(list_xi2[:,N.newaxis]+self.mu-list_eta_hw[N.newaxis,:],self.kernel_Gamma_width)
                

                Fermi2 = function_fermi_occupation(list_xi2[:,N.newaxis]-list_eta_hw[N.newaxis,:],0.,self.beta)
                Fermi1 = function_fermi_occupation(list_xi2[:,N.newaxis],0.,self.beta)
                dFermi = Fermi2-Fermi1  

                for n1, list_epsilon1 in zip([0,1],list_epsilon_k):

                    list_xi1 = N.real(  list_epsilon1 - self.mu )

                    den12  = self.cutoff_denominator(list_xi1[:,N.newaxis] - list_xi2[:,N.newaxis] +list_eta_hw[N.newaxis,:], eta*self.delta_width)


                    for n3, list_epsilon3 in zip([0,1],list_epsilon_kq):

                        list_xi3 = N.real(  list_epsilon3 - self.mu )

                        den32  = self.cutoff_denominator(list_xi3[:,N.newaxis] - list_xi2[:,N.newaxis] +list_eta_hw[N.newaxis,:], eta*self.delta_width)

                        key   = (n1,n2,n3)
                        index = self.index_dictionary[key]

                        #singular_contribution = -eta*dFermi*den12*den32*KR
                        # The previous expression *seems* to have an error, forgetting the contribution from KI. 
                        singular_contribution = -eta*dFermi*den12*den32*(KR+1j*eta*KI)

                        # add the "smooth" part to the I function
                        self.I[index,:,:] += self.conversion_factor*singular_contribution

        return 

    def build_I_singular_intraband_only(self,wedge):
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

            list_xi2 = N.real(  list_epsilon2 - self.mu )

            for eta in [-1.,1.]:

                list_eta_hw =  eta*self.list_hw

                KR = get_KR(list_xi2[:,N.newaxis]+self.mu-list_eta_hw[N.newaxis,:],self.kernel_Gamma_width)

                Fermi2 = function_fermi_occupation(list_xi2[:,N.newaxis]-list_eta_hw[N.newaxis,:],0.,self.beta)
                Fermi1 = function_fermi_occupation(list_xi2[:,N.newaxis],0.,self.beta)
                dFermi = Fermi2-Fermi1  

                for n1, list_epsilon1 in zip([0,1],list_epsilon_k):

                    list_xi1 = N.real(  list_epsilon1 - self.mu )

                    den12  = self.cutoff_denominator(list_xi1[:,N.newaxis] - list_xi2[:,N.newaxis] +list_eta_hw[N.newaxis,:], eta*self.delta_width)


                    for n3, list_epsilon3 in zip([0,1],list_epsilon_kq):

                        list_xi3 = N.real(  list_epsilon3 - self.mu )

                        den32  = self.cutoff_denominator(list_xi3[:,N.newaxis] - list_xi2[:,N.newaxis] +list_eta_hw[N.newaxis,:], eta*self.delta_width)

                        key   = (n1,n2,n3)
                        index = self.index_dictionary[key]

                        singular_contribution = -eta*dFermi*den12*den32*KR

                        if n1 == 0 and n2 == 0 and n3 == 0:
                            singular_contribution_intraband = singular_contribution 
                        else:
                            singular_contribution_intraband = complex(0.,0.)*N.zeros_like(singular_contribution)


                        # add the "smooth" part to the I function
                        self.I[index,:,:] += self.conversion_factor*singular_contribution_intraband 

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

