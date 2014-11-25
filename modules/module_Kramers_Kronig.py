#================================================================================
#
#                  Module Kramers-Kronig
#               ===========================
#
#   This module implements the Kramers Kronig relations, to go from 
#       the 'imaginary' part to the 'real part'.
#================================================================================

from module_Constants import *
import numpy.fft as FFT



class KramersKronig:

    def __init__(self,list_hw):
        """
        Prepare the scheme, from the given frequency grid.
        """

        self.n_hw               = len(list_hw)
        self.n_convolution_size = 2*self.n_hw 

        self.build_convolution_kernel()
        self.build_kernel()
    

    def build_convolution_kernel(self):
        """
        This subroutine computes the array "convolution_kernel" which contains the
        appropriate information to perform the Kramers-Kronig integral
        through a convolution.

        The discretized Kramers-Kronig integral can be expressed as
        R_i = sum_{j} K_{ij} * I_j

        where the kernel K is actually of the form K_{ij} = k[i-j].
        Thus, the sum can be preformed using a convolution, which can be 
        peformed efficiently with FFT!
        """

        one_on_ipi = -1j/N.pi


        # CAUTION! convolution_size should already have been defined!
        self.convolution_kernel = complex(0.,0.)*N.zeros(self.n_convolution_size)

        # Build the linear kernel in direct space
        for l in N.arange(self.n_hw):
            if l == 0:
                kernel = complex(0.,0.)
            elif l == 1:
                kernel = -complex(1.,0.)*2.0*N.log(2.0)
            else:
                kernel =  complex(1.,0.)*             \
                         (  2.0*l *N.log(N.abs(l))    \
                          - (1.+l)*N.log(N.abs(1.+l)) \
                          + (1.-l)*N.log(N.abs(1.-l)))

            self.convolution_kernel[l]     = one_on_ipi * kernel

        for l in N.arange(self.n_hw+1,self.n_convolution_size):
            self.convolution_kernel[l] = -self.convolution_kernel[self.n_convolution_size-l]


    def build_kernel(self):
        self.kernel = FFT.fft(self.convolution_kernel)


    def apply_kramers_kronig_FFT_convolution(self,list_Fi_w):
        """
        This subroutine performs the convolution using FFT.
        """

        list_imaginary = complex(0.,0.)*N.zeros(self.n_convolution_size)

        list_imaginary[:self.n_hw] = list_Fi_w

        KI = FFT.fft(list_imaginary)

        KIF = KI*self.kernel

        list_Fr_w = FFT.ifft(KIF)[:self.n_hw]

        return list_Fr_w 


class KramersKronig_Gamma(KramersKronig):

    def __init__(self,list_hw,Gamma_width):
        """
        Prepare the scheme, from the given frequency grid.
        We assume the width Gamma is a real, positive parameter.
        We also assume the grid is regularly spaced.
        """

        self.n_hw               = len(list_hw)
        self.Delta              = list_hw[2]-list_hw[1]
        self.n_convolution_size = 2*self.n_hw 
        self.gamma              = Gamma_width/self.Delta

        self.generate_kernel()

        def generate_kernel(self):
            self.build_convolution_kernel()
            self.build_kernel()

        def build_convolution_kernel(self):
            """
            This subroutine computes the array "convolution_kernel" which contains the
            appropriate information to perform the Kramers-Kronig integral
            through a convolution.

            The discretized Kramers-Kronig integral can be expressed as
            R_i = sum_{j} K_{ij} * I_j

            where the kernel K is actually of the form K_{ij} = k[i-j].
            Thus, the sum can be preformed using a convolution, which can be 
            peformed efficiently with FFT!
            """

            one_on_ipi = -1j/N.pi


            # Build the linear kernel in direct space

            m_positive = N.arange(0,self.n_hw+1)
            m_negative = N.arange(1-self.n_hw,0)

            m    = N.concatenate([ m_positive , m_negative ])
            mg   = m + 1j*self.gamma
            mg1  = mg+complex(1.,0.)
            mgm1 = mg-complex(1.,0.)

            self.convolution_kernel = one_on_ipi*(mg1*(N.log(-mg)-N.log(-mg1))+mgm1*(N.log(mg)-N.log(mgm1)))



class KramersKronig_KP(KramersKronig_Gamma):
    def generate_kernel(self):
        self.build_convolution_kernel()
        self.convolution_kernel = complex(0.,1.)*N.imag(self.convolution_kernel)
        self.build_kernel()



class KramersKronig_Kd(KramersKronig_Gamma):
    def generate_kernel(self):
        self.build_convolution_kernel()
        self.convolution_kernel = complex(1.,0.)*N.real(self.convolution_kernel)
        self.build_kernel()

