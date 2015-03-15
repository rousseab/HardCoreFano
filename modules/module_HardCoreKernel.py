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


def get_KR(list_energies,kernel_Gamma_width):

    KR = complex(1.,0.)*N.real( get_gamma(list_energies+1j*kernel_Gamma_width)) 

    return KR

def get_KI(list_energies,kernel_Gamma_width):


    limit =  complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff*kernel_Gamma_width/D_cutoff

    KI = -complex(1.,0.)*N.imag( get_gamma(list_energies+1j*kernel_Gamma_width)) - limit

    return KI


class ScatteringKernel:
    """
    The Kramers-Kronig integral will be performed on a dense mesh, and a spline will be used
    to interpolate function to any desired value.

    A finite width will be assumed for the Green function as well as for the scattering Kernel.
    The integrals considered are

    int d xi/(i pi) f(xi) KR(xi+mu)/(xi-(x+i Green_Gamma))
    int d xi/(i pi) f(xi) KI(xi+mu)/(xi-(x+i Green_Gamma))

    int d xi/(i pi) df(xi)/dxi KR(xi+mu)/(xi-(x+i Green_Gamma))
    int d xi/(i pi) df(xi)/dxi KI(xi+mu)/(xi-(x+i Green_Gamma))

    """

    def __init__(self,mu,beta, kernel_Gamma_width,Green_Gamma_width):
        """
        Build the Kramers-Kronig integral
        """

        # chemical potential
        self.mu    = mu

        # inverse temperature
        self.beta  = beta

        # broadening, to avoid singularities in scattering kernel
        self.kernel_Gamma_width = kernel_Gamma_width

        # Green function broadening, to account for lifetime effects
        self.Green_Gamma_width = Green_Gamma_width

        # let's hardcode this parameter for now
        self.spline_order = 3


    def build_scattering_kernel(self):
        """
        Compute the scattering kernel using the KK relation.            
        """

        # Hardcode these parameters for now.

        # build a regular linear energy grid
        n_u   = 4001
        #n_u   = 401
        u_max = 4.*D_cutoff
        u_min =-4.*D_cutoff
        d_u   = (u_max-u_min)/(n_u-1.)
        iloop  = N.arange(n_u)
        list_u = u_min+d_u*iloop

        # Prepare the Kramers-Kronig object
        KK = KramersKronig_Gamma(list_u,self.Green_Gamma_width)

        # Compute the  kernel
        KR = get_KR(list_u+self.mu,self.kernel_Gamma_width)
        KI = get_KI(list_u+self.mu,self.kernel_Gamma_width)

        # Compute the Fermi occupation factor. Note that the chemical potential should NOT be
        # passed to this function 
        f      = function_fermi_occupation(list_u,0.,self.beta)
        df_dxi = d_Fermi_dxi(list_u,0.,self.beta)

        # apply the Kramers-Kronig transformations 
        self.list_KK_fKR  = KK.apply_kramers_kronig_FFT_convolution(f*KR)
        self.list_KK_fKI  = KK.apply_kramers_kronig_FFT_convolution(f*KI)
        self.list_KK_dfKR = KK.apply_kramers_kronig_FFT_convolution(df_dxi*KR)
        self.list_KK_dfKI = KK.apply_kramers_kronig_FFT_convolution(df_dxi*KI)

        self.list_x    = list_u

    def build_and_write_spline_parameters(self,filename):
        """
        Building the parameters for a spline is quite time consuming. It is best to 
        do this once, write to file, and be done with it.

        This routine assumes self.build_scattering_kernel has already created the 
        appropriate arrays.
        """

        # prepare the spline
        self.splmake_tuple_dict = {}

        self.splmake_tuple_dict['Re_fKR']  = splmake(self.list_x, N.real(self.list_KK_fKR), order = self.spline_order)
        self.splmake_tuple_dict['Im_fKR']  = splmake(self.list_x, N.imag(self.list_KK_fKR), order = self.spline_order)

        self.splmake_tuple_dict['Re_dfKR'] = splmake(self.list_x, N.real(self.list_KK_dfKR), order = self.spline_order) 
        self.splmake_tuple_dict['Im_dfKR'] = splmake(self.list_x, N.imag(self.list_KK_dfKR), order = self.spline_order)

        self.splmake_tuple_dict['Re_fKI']  = splmake(self.list_x, N.real(self.list_KK_fKI), order = self.spline_order) 
        self.splmake_tuple_dict['Im_fKI']  = splmake(self.list_x, N.imag(self.list_KK_fKI), order = self.spline_order) 

        self.splmake_tuple_dict['Re_dfKI'] = splmake(self.list_x, N.real(self.list_KK_dfKI), order = self.spline_order) 
        self.splmake_tuple_dict['Im_dfKI'] = splmake(self.list_x, N.imag(self.list_KK_dfKI), order = self.spline_order) 

        write_splmake(self.splmake_tuple_dict, self.spline_order, filename, self.mu, self.beta,\
                                            self.kernel_Gamma_width, self.Green_Gamma_width)

    def read_spline_parameters(self,filename):
        """
        Read in the spline parameters, which must have been generated previously.
        """

        # prepare the spline
        self.splmake_tuple_dict, file_mu, file_beta, \
        file_kernel_Gamma_width, file_Green_Gamma_width  = read_splmake(filename)

        tol = 1e-8
        if N.abs(file_mu   - self.mu) > tol or \
            N.abs(file_beta - self.beta) > tol or \
              N.abs(file_Green_Gamma_width - self.Green_Gamma_width) > tol or \
                N.abs(file_kernel_Gamma_width - self.kernel_Gamma_width) > tol:
            print 'ERROR: the splmake_tuple file contains parameters inconsistent with this calculation'
            print 'EXITING NOW'
            sys.exit()

        return

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
