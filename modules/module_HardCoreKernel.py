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
    This class will essentially produce the scattering kernel function L, which is given by

    L(xi,nhw) = f(xi) KR(xi+mu-nhw) + i KramerKronig (KI)(xi+mu-nhw)

    The Kramers-Kronig integral will be performed on a dense mesh, and a spline will be used
    to interpolate function to any desired value.
    """

    def __init__(self,mu,beta, kernel_Gamma_width):
        """
        Build the Kramers-Kronig integral
        """

        # chemical potential
        self.mu    = mu

        # inverse temperature
        self.beta  = beta

        # broadening, to avoid singularities in scattering kernel
        self.kernel_Gamma_width = kernel_Gamma_width

    def build_scattering_kernel(self):
        """
        Compute the scattering kernel using the KK relation.            
        """

        # Hardcode these parameters for now.

        # build a regular linear energy grid
        n_u   = 4001
        u_max = 4.*D_cutoff
        u_min =-4.*D_cutoff
        d_u   = (u_max-u_min)/(n_u-1.)
        iloop  = N.arange(n_u)
        list_u = u_min+d_u*iloop

        # Prepare the Kramers-Kronig object
        KK = KramersKronig(list_u)

        # Compute the anti-symmetric kernel
        KI = get_KI(list_u+self.mu,self.kernel_Gamma_width)

        # Compute the Fermi occupation factor. Note that the chemical potential should NOT be
        # passed to this function 
        f  = function_fermi_occupation(list_u,0.,self.beta)

        # the integrand in the Kramers-Kronig integral
        list_fKI = f*KI

        # apply the Kramers-Kronig transformation to the anti-symmetric kernel
        self.list_KK_fKI = N.real(1j*KK.apply_kramers_kronig_FFT_convolution(list_fKI))
        self.list_x    = list_u

    def build_and_write_spline_parameters(self,filename):
        """
        Building the parameters for a spline is quite time consuming. It is best to 
        do this once, write to file, and be done with it.

        This routine assumes self.build_scattering_kernel has already created the 
        appropriate arrays.
        """

        # prepare the splne
        self.splmake_tuple = splmake(self.list_x, self.list_KK_fKI )
        write_splmake(self.splmake_tuple,filename,self.mu,self.beta,self.kernel_Gamma_width)

    def read_spline_parameters(self,filename):
        """
        Read in the spline parameters, which must have been generated previously.
        """

        # prepare the splne
        self.splmake_tuple, file_mu, file_beta, file_kernel_Gamma_width = read_splmake(filename)

        tol = 1e-8
        if N.abs(file_mu   - self.mu) > tol or \
            N.abs(file_beta - self.beta) > tol or \
              N.abs(file_kernel_Gamma_width - self.kernel_Gamma_width) > tol:
            print 'ERROR: the splmake_tuple file contains parameters inconsistent with this calculation'
            print 'EXITING NOW'
            sys.exit()

        return

    def get_LR_oneD(self,list_xi):
        """
        Compute the contribution coming from the real kernel
        """

        f = function_fermi_occupation(list_xi,0.,self.beta)

        # Compute the symmetric kernel
        KR = get_KR(list_xi,self.kernel_Gamma_width)

        LR = complex(1.,0.)*f*KR    

        return LR

    def get_LR_twoD(self,list_xi,list_eta_hw):
        """
        Compute the contribution coming from the real kernel
        """

        fermi_dirac  = function_fermi_occupation(list_xi,0.,self.beta)

        u  = list_xi[:,N.newaxis] - list_eta_hw[N.newaxis,:] + self.mu
        f  = fermi_dirac[:,N.newaxis]


        KR = get_KR(u,self.kernel_Gamma_width)

        LR = complex(1.,0.)*f*KR

        return LR

    def get_LI_oneD(self,list_xi):
        """
        Compute the contribution coming from the imaginary kernel
        """

        # evaluate spline            
        u = N.real ( list_xi )

        LI = complex(1.,0.)*spleval(self.splmake_tuple, u)

        return LI

    def get_LI_twoD(self,list_xi,list_eta_hw):
        """
        Compute the contribution coming from the imaginary kernel
        """
        # evaluate spline            
        u = N.real ( list_xi[:,N.newaxis]  - list_eta_hw[N.newaxis,:] )

        LI = complex(1.,0.)*spleval(self.splmake_tuple, u)

        return LI

    def get_L_oneD(self,list_xi):
        """
        Compute the full contribution to the scattering kernel
        """

        L = self.get_LR_oneD(list_xi)+self.get_LI_oneD(list_xi)

        return L

    def get_L_twoD(self,list_xi,list_eta_hw):
        """
        Compute the full contribution to the scattering kernel
        """

        L = self.get_LR_twoD(list_xi,list_eta_hw)+self.get_LI_twoD(list_xi,list_eta_hw)

        return L

