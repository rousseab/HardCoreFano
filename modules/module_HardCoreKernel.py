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

gamma_constant =  -2.*N.pi*hvF**2/D_cutoff

def get_gamma(list_epsilon):
    """
    The argument should be epsilon_{nk} = xi_{nk}+mu, in eV.

    The returned value is in eV x a0^2
    """
    gamma = complex(1.,0.)*gamma_constant*list_epsilon/D_cutoff

    return gamma

def get_gammaA(list_z):
    """
    The "real" kernel function is given by

        gamma(z) =  2 pi (hvF)^2/D ( D/z 1/ln[1-D^2/z^2])

    This function is wildly non-analytic near zero. We thus wish to
    extract its behavior at inifnity.    

    z-> infinity limit:
    gamma_inf (z) =- 2 pi (hvF)^2/D  z/D

    Pade approximation of the leftover, to first order:
    gamma_1   (z) =  2 pi (hvF)^2/D  ( 3x^3+2x)/(6x^4+3x^2+2),  x  = z/D

    gammaA = gamma_inf + gamma_1 

           = 2 pi (hvF)^2/D [-  x^5/(x^4 + 1/2 x^2 + 1/3)]

    The argument of the function should always be offset by mu!

    The returned value is in eV x a0^2
    """

    prefactor =  complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff

    x = list_z/D_cutoff
    x2= x*x
    x3= x2*x
    x4= x3*x
    x5= x4*x


    numerator   = -x5
    denominator = complex(1.,0.)*x4+complex(0.5,0.)*x2+complex(1.,0.)/complex(3.,0)

    gamma_A   = prefactor*numerator/denominator 

    return gammaA


def get_gammaA_POLES_and_RESIDUES():
    """

    gammaA = 2 pi (hvF)^2/D [-  x^5/(x^4 + 1/2 x^2 + 1/3)]

    The poles are given by  Z_j = rho e^{i theta_j},  rho = (1/3)^{1/4} D,

        theta_j = { pi/2-phi, pi/2+phi, 3pi/2 - phi, 3pi/2 + phi}, phi = 1/2 arctan(sqrt{13/3})

    """

    prefactor =  complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff

    phi = 0.5*N.arctan(N.sqrt(13./3.))

    rho = (1./3.)**(0.25)

    list_angles =  N.array([ N.pi/2-phi, N.pi/2+phi, 3*N.pi/2-phi, 3*N.pi/2+phi])
    list_phases =  N.exp(1j*list_angles) 

    list_x  =  rho*list_phases 

    list_POLES  =  D_cutoff*list_x

    list_RESIDUES = []

    for i, x_i in enumerate(list_x):

        numerator = -x_i**5

        denominator = complex(1.,0.)

        for j, x_j in enumerate(list_x):
            if j != i:
                denominator = denominator*(x_i-x_j)

        residue = prefactor*numerator/denominator

        list_RESIDUES.append(residue)

    list_RESIDUES = N.array(list_RESIDUES)

    return list_POLES, list_RESIDUES 

def get_gammaC(list_z):
    """
    The "real" kernel function is given by

        gamma(z) =  2 pi (hvF)^2/D ( D/z 1/ln[1-D^2/z^2])

    This function is wildly non-analytic near zero. 

    gammaC = gamma - gammaA

    The argument of the function should always be offset by mu!

    The returned value is in eV x a0^2
    """

    prefactor =  complex(1.,0.)*2.*N.pi*hvF**2/D_cutoff

    gammaA = get_gammaA(list_z)

    x = complex(1.,0.)*list_z/D_cutoff
    x2= x*x

    denominator = x*N.log(complex(1.,0.)-complex(1.,0.)/x2)

    gamma = prefactor/denominator 

    gammaC = gamma - gammaA   

    return gammaC


class ScatteringKernel:
    """
    This class will essentially produce the scattering kernel function L, which
    can be decomposed as

    L(xi,hw)  = LA1(xi,hw) + LA2(xi,hw) + LC(xi,hw)

    LA1(xi,hw) = sum_{i} R_i f(z_i-mu)/(z_i-mu-(xi-hw)),  {R_i,z_i} residues/poles of gamma_A

    LA2(xi,hw) = f(xi) gammaA(xi+mu-hw) # Careful! hw does not appear in f()

    LC(xi,hw) = f(xi-hw) KCR(xi+mu-hw) + i KramerKronig (KCI)(xi+mu-hw)

    The Kramers-Kronig integral will be performed on a dense mesh, and a spline will be used
    to interpolate function to any desired value.
    """

    def __init__(self,mu,beta, delta):
        """
        Build the Kramers-Kronig integral
        """

        # chemical potential
        self.mu    = mu

        # inverse temperature
        self.beta  = beta

        # broadening, to avoid singularities in scattering kernel
        self.delta = delta

        # save the poles and residues gammaA
        self.list_POLES, self.list_RESIDUES = get_gammaA_POLES_and_RESIDUES()

        # this is the combination that will be used
        self.ZA_mu = self.list_POLES-self.mu

        # complex Fermi-Dirac
        f = complex(1.,0.)/ ( N.exp(self.beta*self.ZA_mu)+complex(1.,0.) )

        self.RAf = self.list_RESIDUES*f

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

        # Compute symmetric and anti-symmetric kernels
        gC_plus    = get_gammaC(list_u+self.mu+1j*self.delta)
        gC_minus   = get_gammaC(list_u+self.mu-1j*self.delta)

        KCI = (gC_minus-gC_plus)/(2.*1j)
        KCR = (gC_minus+gC_plus)/2.

        # Compute the Fermi occupation factor. Note that the chemical potential should NOT be
        # passed to this function 
        f  = function_fermi_occupation(list_u,0.,self.beta)

        # the symmetric part of the scattering kernel
        #  Only take real part (imaginary part should be zero anyways) because
        # the spline interpolation requires floats.
        list_fKCR = N.real(f*KCR)

        # the integrand in the Kramers-Kronig integral
        list_fKCI = f*KCI

        # apply the Kramers-Kronig transformation to the anti-symmetric kernel
        list_KK_fKCI = N.real(1j*KK.apply_kramers_kronig_FFT_convolution(list_fKCI))

        # sum the two, as it appears some of the singularities present in each term ends up smoothed out in the sum
        self.LC_kernel = list_fKCR + list_KK_fKCI
        self.list_x    = list_u

    def build_and_write_spline_parameters(self,filename):
        """
        Building the parameters for a spline is quite time consuming. It is best to 
        do this once, write to file, and be done with it.

        This routine assumes self.build_scattering_kernel has already created the 
        appropriate arrays.
        """

        # prepare the splne
        self.splmake_tuple = splmake(self.list_x, self.LC_kernel )
        write_splmake(self.splmake_tuple,filename,self.mu,self.beta,self.delta)

    def read_spline_parameters(self,filename):
        """
        Read in the spline parameters, which must have been generated previously.
        """

        # prepare the splne
        self.splmake_tuple, file_mu, file_beta, file_delta = read_splmake(filename)

        tol = 1e-8
        if N.abs(file_mu   - self.mu) > tol or \
            N.abs(file_beta - self.beta) > tol or \
              N.abs(file_delta- self.delta) > tol:
            print 'ERROR: the splmake_tuple file contains parameters inconsistent with this calculation'
            print 'EXITING NOW'
            sys.exit()

        return

    def get_LA1(self,list_xi,list_hw=None):
        """
        Compute the LA1 contribution
        """

        if list_hw == None:
            u = list_xi
        else:            
            u = list_xi[:,N.newaxis]  - list_hw[N.newaxis,:]

        LA1 =   self.RAf[0]/(self.ZA_mu[0]-u) +\
                self.RAf[1]/(self.ZA_mu[1]-u) +\
                self.RAf[2]/(self.ZA_mu[2]-u) +\
                self.RAf[3]/(self.ZA_mu[3]-u)


        return LA1

    def get_LA2(self,list_xi,list_hw=None):
        """
        Compute the LA2 contribution
        """

        fermi_dirac  = function_fermi_occupation(list_xi,0.,self.beta)

        if list_hw == None:
            u  = list_xi
            f  = fermi_dirac
        else:
            u  = list_xi[:,N.newaxis] - list_hw[N.newaxis,:]
            f  = fermi_dirac[:,N.newaxis]

        gammaA = get_gammaA(u+self.mu )

        LA2 = f*gammaA

        return LA2

    def get_LC(self,list_xi,list_hw=None):
        """
        Compute the gammaC contribution to the scattering kernel
        """

        # evaluate spline            
        if list_hw == None:
            u = N.real ( list_xi )
        else:            
            u = N.real ( list_xi[:,N.newaxis]  - list_hw[N.newaxis,:] )

        LC = complex(1.,0.)*spleval(self.splmake_tuple, u)

        return LC

    def get_L(self,list_xi,list_hw=None):
        """
        Compute the full contribution to the scattering kernel
        """

        if list_hw == None:
            L = self.get_LC(list_xi)+self.get_LA1(list_xi)+self.get_LA2(list_xi)
        else:
            L = self.get_LC(list_xi,list_hw)+self.get_LA1(list_xi,list_hw)+self.get_LA2(list_xi,list_hw)

        return L

