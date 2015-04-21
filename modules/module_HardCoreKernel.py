#================================================================================
#
#           Module HardCoreKernel
#           =====================
#
#   This module implements the simplest analytical functions which describe
#   the hard core kernel for graphene.
#
#================================================================================
import sys
from module_Constants import *
from module_Functions import *
from module_NETCDF    import *

from module_Kramers_Kronig import *

from scipy.interpolate import splmake, spleval, spline


# Define a few useful constants which are only useful for the Hard Core Kernel

k_cutoff = N.sqrt(2.*N.pi/Area) # in 1/bohr
D_cutoff = hvF*k_cutoff         # in eV


def get_gamma(list_z):
    """

    The "real" kernel function is given by

        gamma(z) =  2 pi (hvF)^2/D ( D/z 1/ln[1-D^2/z^2])

    This function is wildly non-analytic near zero. 

    The argument of the function should always be offset by mu!

    The returned value is in eV x a0^2
    """

    prefactor =  complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff

    x = complex(1.,0.)*list_z/D_cutoff
    x2= x*x

    denominator = x*N.log(complex(1.,0.)-complex(1.,0.)/x2)

    gamma = prefactor/denominator 

    return gamma

def get_KR_inf(list_energies):
    KR_inf  =  -complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff**2*list_energies

    return KR_inf

def get_KR(list_energies,kernel_Gamma_width):

    KR = complex(1.,0.)*N.real( get_gamma(list_energies+1j*kernel_Gamma_width)) 

    return KR

def get_SR_inf(list_epsilon,mu,Green_Gamma_width):
    """
    This function returns the effective contribution coming from the 
    linear term in the KR kernel. The integral over this linear part is 
    formally divergent, but the divergence can be regularized. We are
    left with the following contribution.

    See document "computing_J.pdf" for details.
    """

    prefactor =  -complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff**2
    
    Term_1 = N.log(1.+list_epsilon**2/Green_Gamma_width**2)

    Term_2 = 1j*N.pi-2.*1j*N.arctan(list_epsilon/Green_Gamma_width)

    energy_factor =  list_epsilon +mu + 1j*Green_Gamma_width

    SR_inf  =  prefactor/(2.*N.pi*1j)* energy_factor *(Term_1+ Term_2)

    return SR_inf

def get_SR_finite_cutoff_correction(list_epsilon,mu,Lambda_cutoff,Green_Gamma_width):
    """
    The integral

    SR_min_inf(epsilon,Gamma) =  int dxi/(i pi) f(xi) K^R(xi+mu) - K^{infty}(xi+mu)/(xi-(epsilon + i Gamma) )

    Converges very slowly with the integration cutoff Lambda_cutoff. Since Kramers-Kronig necessarily
    uses such a cutoff, this routine implements the first correction to the cutoff integral when the
    integration bounds go to infinity.

    
    See document "computing_J.pdf" for details.
    """
    prefactor = -1j*hvF**2/2.        
    den = 1./(list_epsilon+mu+1j*Green_Gamma_width)

    Term_1 = N.log( (Green_Gamma_width**2+(Lambda_cutoff+list_epsilon)**2)/(Lambda_cutoff-mu)**2)
    Term_2 = 1j*N.pi-1j*2*N.arctan((Lambda_cutoff+list_epsilon)/Green_Gamma_width)

    correction = prefactor*den*(Term_1+Term_2)

    return correction 


def get_KI(list_energies,kernel_Gamma_width):


    limit =  complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff*kernel_Gamma_width/D_cutoff

    KI = -complex(1.,0.)*N.imag( get_gamma(list_energies+1j*kernel_Gamma_width)) - limit

    return KI

class ScatteringKernel:
    """
    This oject computes and stores  generalized Kramers-Kronig integrals, as well as the
    the necessary technology to interpolate to arbitrary values of the arguments.

    A finite width will be assumed for the Green function as well as for the scattering Kernel.

    The integrals considered are

    SR(x, kappa) = int d xi/(i pi) f(xi) KR(xi+mu)/(xi-(x+i kappa Green_Gamma))
    SI(x, kappa) = int d xi/(i pi) f(xi) KI(xi+mu)/(xi-(x+i kappa Green_Gamma))

    TR(x, kappa, eta hw) = int d xi/(i pi) ( f(xi) - f(xi+eta hw) ) KR(xi+mu)/(xi-(x+i kappa Green_Gamma))
    TI(x, kappa, eta hw) = int d xi/(i pi) ( f(xi) - f(xi+eta hw) ) KI(xi+mu)/(xi-(x+i kappa Green_Gamma))

    Note that SR is formally divergent, so various regularization schemes are introduced.

    """

    def __init__(self,mu,beta, kernel_Gamma_width, Green_Gamma_width, spline_order = 3):
        """
        Prepare the basics for this object. It can be invoked to read or write the
        integrals.        
        """

        # chemical potential
        self.mu    = mu

        # inverse temperature
        self.beta  = beta

        # broadening, to avoid singularities in scattering kernel
        self.kernel_Gamma_width = kernel_Gamma_width

        # Green function broadening, to account for lifetime effects
        self.Green_Gamma_width = Green_Gamma_width

        self.spline_order = spline_order 

    def build_scattering_kernel_integrals(self, list_hw_ext, n_xi_grid = 4001):
        """
        Compute the various integrals
        """
        # External frequencies at which TA must be computed
        self.list_hw_ext =  list_hw_ext 
        self.list_eta    =  N.array([-1,1])

        # build a regular linear energy grid
        self.n_xi_grid   =  n_xi_grid 

        xi_max = 4.*D_cutoff
        xi_min =-4.*D_cutoff

        Lambda_cutoff = xi_max 

        d_xi   = (xi_max-xi_min)/(self.n_xi_grid-1.)
        iloop  = N.arange(self.n_xi_grid)

        self.list_xi = xi_min+d_xi*iloop

        # Prepare the Kramers-Kronig object
        KK = KramersKronig_Gamma(self.list_xi,self.Green_Gamma_width)

        # Compute the  kernels
        KR = get_KR(self.list_xi+self.mu,self.kernel_Gamma_width)
        KI = get_KI(self.list_xi+self.mu,self.kernel_Gamma_width)

        KR_inf = get_KR_inf(self.list_xi+self.mu)

        # Compute the Fermi occupation factor. Note that the chemical potential should NOT be
        # passed to this function 
        f = function_fermi_occupation(self.list_xi,0.,self.beta)

        # Compute the Heaviside plateau function.
        Theta = Heaviside_plateau(self.list_xi)


        # apply the Kramers-Kronig transformations 
        print ' ====   Computing SR and SI ====='
        self.SI  = KK.apply_kramers_kronig_FFT_convolution(f*KI)

        # Compute SR in multiple steps, to avoid singular kernel

        
        # Keep the following arrays in memory, for testing and debugging!            
        self.SR_1_cutoff = KK.apply_kramers_kronig_FFT_convolution(f*(KR-KR_inf)) # converges slowly because of cutoff
        self.SR_1_cutoff_correction = \
                get_SR_finite_cutoff_correction(self.list_xi,self.mu,Lambda_cutoff,self.Green_Gamma_width)

        self.SR_1 = self.SR_1_cutoff+ self.SR_1_cutoff_correction 

        self.SR_2 = KK.apply_kramers_kronig_FFT_convolution((f-Theta)*KR_inf)
        self.SR_inf = get_SR_inf(self.list_xi,self.mu,self.Green_Gamma_width)

        self.SR  = self.SR_1+self.SR_2+self.SR_inf 

        # Initialize arrays with dimensions [ eta, hw, n_xi_grid]
        # Arrays are C ordered by default, so last index loops fastest.
        self.TR = complex(0.,0.)*N.zeros([2,len(self.list_hw_ext),self.n_xi_grid])
        self.TI = complex(0.,0.)*N.zeros([2,len(self.list_hw_ext),self.n_xi_grid])

        # Bite the bullet, do this fairly heavy computation!
        for i_eta, eta in enumerate(self.list_eta):
            for i_hw_ext, hw_ext in enumerate( self.list_hw_ext ):
        
                print ' ====   Computing TR and TI, eta = %i,  hw_ext = %4.3f meV ====='%(i_eta,1000*hw_ext)
                # Compute the Fermi occupation factor, shifted by frequency
                # 
                f_ehw = function_fermi_occupation(self.list_xi+eta*hw_ext,0.,self.beta)

                df = f-f_ehw

                self.TR[i_eta,i_hw_ext,:] = KK.apply_kramers_kronig_FFT_convolution(df*KR)
                self.TI[i_eta,i_hw_ext,:] = KK.apply_kramers_kronig_FFT_convolution(df*KI)

        return

    def build_and_write_spline_parameters(self,filename):
        """
        Building the parameters for a spline is quite time consuming. It is best to 
        do this once, write to file, and be done with it.

        This routine assumes self.build_scattering_kernel has already created the 
        appropriate arrays.
        """

        # generate the splines

        # The function splmake has inteface:
        #
        #           (xk, Ck, order)  = splmake(xk, yk, order)   
        #
        # There is no point in storing "xk" and "order" as part of the output tuple. 
        # Spleval will need this information in the form of a tuple, which we can re-build
        #  at evalation time.

        print ' ====   Computing spline of SR ====='
        list_xi, re_SR_spline, order = splmake(self.list_xi, N.real(self.SR), order = self.spline_order)
        list_xi, im_SR_spline, order = splmake(self.list_xi, N.imag(self.SR), order = self.spline_order)

        print ' ====   Computing spline of SI ====='
        list_xi, re_SI_spline, order = splmake(self.list_xi, N.real(self.SI), order = self.spline_order)
        list_xi, im_SI_spline, order = splmake(self.list_xi, N.imag(self.SI), order = self.spline_order)

        # same length for all splines
        spline_length = len(im_SI_spline)


        re_TR_spline = N.zeros([2,len(self.list_hw_ext),spline_length ])
        im_TR_spline = N.zeros([2,len(self.list_hw_ext),spline_length ])
        re_TI_spline = N.zeros([2,len(self.list_hw_ext),spline_length ])
        im_TI_spline = N.zeros([2,len(self.list_hw_ext),spline_length ])

        for i_eta, eta in enumerate(self.list_eta):
            for i_hw_ext, hw_ext in enumerate( self.list_hw_ext):
                print ' ====   Computing spline for TR and TI, eta = %i,  hw_ext = %4.3f meV ====='%(eta,1000*hw_ext)

                list_xi, re_TR_spline[i_eta,i_hw_ext,:], order = splmake(self.list_xi, N.real(self.TR[i_eta,i_hw_ext,:]), order = self.spline_order)
                list_xi, im_TR_spline[i_eta,i_hw_ext,:], order = splmake(self.list_xi, N.imag(self.TR[i_eta,i_hw_ext,:]), order = self.spline_order)

                list_xi, re_TI_spline[i_eta,i_hw_ext,:], order = splmake(self.list_xi, N.real(self.TI[i_eta,i_hw_ext,:]), order = self.spline_order)
                list_xi, im_TI_spline[i_eta,i_hw_ext,:], order = splmake(self.list_xi, N.imag(self.TI[i_eta,i_hw_ext,:]), order = self.spline_order)

        write_splmake(  re_SR_spline, im_SR_spline, re_SI_spline, im_SI_spline,
                        re_TR_spline, im_TR_spline, re_TI_spline, im_TI_spline,
                        filename, self.spline_order, self.list_xi, self.list_hw_ext, 
                        self.mu, self.beta, self.kernel_Gamma_width, self.Green_Gamma_width)

    def read_spline_parameters(self,filename):
        """
        Read in the spline parameters, which must have been generated previously.
        """

        # prepare the spline
        self.re_SR_spline, self.im_SR_spline, self.re_SI_spline, self.im_SI_spline,\
        self.re_TR_spline, self.im_TR_spline, self.re_TI_spline, self.im_TI_spline,\
        self.list_xi, self.list_hw_ext, file_mu, file_beta, \
        file_kernel_Gamma_width, file_Green_Gamma_width, file_spline_order = read_splmake(filename)

        tol = 1e-8

        test_mu    = N.abs(file_mu   - self.mu) > tol
        test_beta  = N.abs(file_beta - self.beta) > tol 
        test_Green = N.abs(file_Green_Gamma_width - self.Green_Gamma_width)> tol
        test_kernel= N.abs(file_kernel_Gamma_width - self.kernel_Gamma_width) > tol
        test_order = N.abs(file_spline_order - self.spline_order) > tol

        if  test_mu or test_beta or test_Green or test_kernel or test_order: 
            print 'ERROR: the splmake file contains parameters inconsistent with this calculation'
            print 'EXITING NOW'
            sys.exit()

        return

    def get_spline_SR(self, list_epsilon, sign_Gamma):
        """
        Compute the interpolated value for SR, using splines, which we assume are already computed.
        """
        # evaluate spline            
        u = N.real ( list_epsilon )

        re_tuple = (self.list_xi, self.re_SR_spline, self.spline_order)
        im_tuple = (self.list_xi, self.im_SR_spline, self.spline_order)

        S_R = sign_Gamma*complex(1.,0.)*spleval(re_tuple,u)+ complex(0.,1.)*spleval(im_tuple,u)

        return S_R

    def get_spline_SI(self, list_epsilon, sign_Gamma):
        """
        Compute the interpolated value for SR, using splines, which we assume are already computed.
        """
        # evaluate spline            
        u = N.real ( list_epsilon )

        re_tuple = (self.list_xi, self.re_SI_spline, self.spline_order)
        im_tuple = (self.list_xi, self.im_SI_spline, self.spline_order)

        S_I = sign_Gamma*complex(1.,0.)*spleval(re_tuple,u)+ complex(0.,1.)*spleval(im_tuple,u)

        return S_I



    def get_integral_KK(self,term_name, list_xi, sign_Gamma):
        """
        Compute the contribution coming from the imaginary kernel.
        """

        # evaluate spline            
        u = N.real ( list_xi )

        f_KK = \
                sign_Gamma*complex(1.,0.)*spleval(self.splmake_tuple_dict['Re_%s'%term_name],u)+\
                           complex(0.,1.)*spleval(self.splmake_tuple_dict['Im_%s'%term_name],u)

        return f_KK 


