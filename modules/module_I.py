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

        self.mu  = mu
        self.beta= beta

        # Create the scattering kernel object
        self.SK = ScatteringKernel(self.mu,self.beta,self.delta_width)

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

        Y(xi1,xi2,eta_hw) =  ( L(xi1,0) -L(xi2-eta_hw, 0 ) ) / (xi1-[xi2-eta_hw]).

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
        tol = 1e-5

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

        list_KR = []
        for eta in [-1.,1.]:

            list_eta_hw =  eta*N.real(self.z)

            list_KR.append( N.real( get_gamma(self.mu-list_eta_hw+1j*self.delta_width) ) )



        for n2, list_epsilon2 in zip([0,1],list_epsilon_k):

            list_xi2 = N.real(  list_epsilon2 - self.mu )

            delta_xi2 =  -N.imag ( self.cutoff_denominator(list_xi2,fadeeva_width=self.delta_width) )/N.pi

            for eta, KR in zip([-1.,1.],list_KR):

                list_eta_hw =  eta*N.real(self.z)

                WKR = N.real(self.z)*KR

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

                    if eta > 0.:
                        den1 = self.cutoff_denominator(list_xi1[:,N.newaxis]+list_eta_hw[N.newaxis,:],fadeeva_width=self.delta_width)
                    else:
                        den1 = N.conjugate( self.cutoff_denominator(list_xi1[:,N.newaxis]+list_eta_hw[N.newaxis,:],fadeeva_width=self.delta_width) )

                    for n3, list_epsilon3, Y32 in zip([0,1],list_epsilon_kq,list_Y32):

                        list_xi3 = N.real(  list_epsilon3 - self.mu )

                        if eta > 0.:
                            den3 = self.cutoff_denominator(list_xi3[:,N.newaxis]+list_eta_hw[N.newaxis,:],fadeeva_width=self.delta_width)
                        else:
                            den3 = N.conjugate( self.cutoff_denominator(list_xi3[:,N.newaxis]+list_eta_hw[N.newaxis,:],fadeeva_width=self.delta_width) )


                        key   = (n1,n2,n3)
                        index = self.index_dictionary[key]

                        # cutoff the 1/0 
                        den13 = N.real ( self.cutoff_denominator(list_xi1-list_xi3,fadeeva_width=1e-6) )

                        smooth_contribution   = eta*den13[:,N.newaxis]*(Y12-Y32) 

                        singular_contribution = -delta_xi2[:,N.newaxis]*den1*den3*WKR[N.newaxis,:]

                        # add the "smooth" part to the I function
                        self.I[index,:,:] += self.conversion_factor*(smooth_contribution+singular_contribution)

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
        G = complex(1.,0.)/(self.iwm[N.newaxis,:] - list_xi[:,N.newaxis])

        return G

    def get_frequency_double_green_function(self, list_xi):
        """
        Produce the tensor
        D [ nk, hw, iwm]
        D = (hw+idelta)/( (iwm-xik)**2 - (hw+idelta)^2 ) 
        """

        t1 = (self.iwm[N.newaxis,N.newaxis,:]-list_xi[:,N.newaxis,N.newaxis])**2
        t2 = self.z2[N.newaxis,:,N.newaxis]

        D  = self.z[N.newaxis,:,N.newaxis]/(t1-t2)

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

