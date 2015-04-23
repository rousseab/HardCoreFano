#================================================================================
#
# This test module implements necessary tools to compute the J term
# exactly using numerical integration.
#
#================================================================================
import common
reload(common)
from common import *

import scipy.integrate as SI
from functools import partial
# import the usual modules

import numpy as N

class NetcdfHandler(object):

    def __init__(self, filename):
        self.filename = filename
    
    def open_ncfile(self,mode):
        self.ncfile = Dataset(self.filename,mode)

    def close_ncfile(self):
        self.ncfile.close()

    def write_attributes(self,mu, beta, Green_Gamma_width, kernel_Gamma_width, hw_ph, list_xi_1, list_xi_3):

        # --- set various attributes, identifying the parameters of the computation ----
        setattr(self.ncfile,'mu',mu) 
        setattr(self.ncfile,'beta',beta) 
        setattr(self.ncfile,'Green_Gamma_width',Green_Gamma_width) 
        setattr(self.ncfile,'kernel_Gamma_width',kernel_Gamma_width) 
        setattr(self.ncfile,'phonon_frequency',hw_ph) 

        # --- Create dimensions ----
        self.ncfile.createDimension("xi_1",len(list_xi_1))
        self.ncfile.createDimension("xi_3",len(list_xi_3))
     
        # --- Write data ----
        Xi_n1k  = self.ncfile.createVariable("list_xi_n1k",'d',('xi_1',))
        Xi_n3kq = self.ncfile.createVariable("list_xi_n3kq",'d',('xi_3',))

        Xi_n1k[:]  = list_xi_1
        Xi_n3kq[:] = list_xi_3

        return

    def read_attributes(self):

        # --- set various attributes, identifying the parameters of the computation ----
        mu = self.ncfile.mu[0]
        Green_Gamma_width = self.ncfile.Green_Gamma_width[0]
        kernel_Gamma_width = self.ncfile.kernel_Gamma_width[0]
        hw_ph = self.ncfile.phonon_frequency[0]
        beta  = self.ncfile.beta[0]

        list_xi_1 = self.ncfile.variables["list_xi_n1k"][:]
        list_xi_3 = self.ncfile.variables["list_xi_n3kq"][:]

        return list_xi_1, list_xi_3, mu, Green_Gamma_width, kernel_Gamma_width, hw_ph, beta

    def write_J_array(self,J_array,name):

        Re_J = self.ncfile.createVariable("Re_%s"%name,'d',('xi_1','xi_3'))
        Im_J = self.ncfile.createVariable("Im_%s"%name,'d',('xi_1','xi_3'))

        Re_J[:,:] = N.real(J_array)
        Im_J[:,:] = N.imag(J_array)

    def read_J_array(self,name):

        Re_J = self.ncfile.variables["Re_%s"%name][:]
        Im_J = self.ncfile.variables["Im_%s"%name][:]
        J    = complex(1.,0.)*Re_J + complex(0.,1.)*Im_J

        return J

def get_exact_J_integrand(xi, xi_n1k, xi_n2k, xi_n3kq, hw, mu, beta, kernel_Gamma_width, Green_Gamma_width):
    """
    This function computes the exact integrand leading to the J function.
    """

    J_integrand_df = complex(0.,0.)
    J_integrand_f  = complex(0.,0.)

    list_xi = N.array([xi])

    KR = get_KR(list_xi+mu,kernel_Gamma_width)[0]
    KI = get_KI(list_xi+mu,kernel_Gamma_width)[0]

    f_xi = function_fermi_occupation(list_xi,0.,beta)[0]

    for eta in [-1,1]:
        prefactor = eta/(2.*N.pi*1j)

        kernel = KR-1j*eta*KI

        D_n1k  = complex(1.,0.)/(xi - xi_n1k  + 1j*eta*Green_Gamma_width)
        D_n3kq = complex(1.,0.)/(xi - xi_n3kq + 1j*eta*Green_Gamma_width)

        D_n2k_1 = complex(1.,0.)/(xi - xi_n2k -eta*hw - 1j*Green_Gamma_width)
        D_n2k_2 = complex(1.,0.)/(xi - xi_n2k -eta*hw + 1j*Green_Gamma_width)

        D_n2k_3 = complex(1.,0.)/(xi - xi_n2k - hw + 1j*eta*Green_Gamma_width)
        D_n2k_4 = complex(1.,0.)/(xi - xi_n2k + hw + 1j*eta*Green_Gamma_width)

        f_xi_ehw = function_fermi_occupation(list_xi-eta*hw,0.,beta)[0]

        Term_1  =  ( f_xi - f_xi_ehw ) * (D_n2k_1 - D_n2k_2) 

        Term_2  =   f_xi * (D_n2k_3 - D_n2k_4) 

        common_factor = prefactor*kernel*D_n1k*D_n3kq 

        J_integrand_df += common_factor*Term_1  
        J_integrand_f  += common_factor*Term_2  

    J_integrand = J_integrand_df+ J_integrand_f  

    return J_integrand, J_integrand_df, J_integrand_f  

def get_exact_linear_kernel_J_integrand(xi, xi_n1k, xi_n2k, xi_n3kq, hw, mu, beta, kernel_Gamma_width, Green_Gamma_width):
    """
    This function computes the exact integrand leading to the J function.
    """

    J_integrand_df = complex(0.,0.)
    J_integrand_f  = complex(0.,0.)

    list_xi = N.array([xi])

    KR_inf = get_KR_inf(list_xi+mu)[0]

    f_xi = function_fermi_occupation(list_xi,0.,beta)[0]

    for eta in [-1,1]:
        prefactor = eta/(2.*N.pi*1j)

        kernel = KR_inf

        D_n1k  = complex(1.,0.)/(xi - xi_n1k  + 1j*eta*Green_Gamma_width)
        D_n3kq = complex(1.,0.)/(xi - xi_n3kq + 1j*eta*Green_Gamma_width)

        D_n2k_1 = complex(1.,0.)/(xi - xi_n2k -eta*hw - 1j*Green_Gamma_width)
        D_n2k_2 = complex(1.,0.)/(xi - xi_n2k -eta*hw + 1j*Green_Gamma_width)

        D_n2k_3 = complex(1.,0.)/(xi - xi_n2k - hw + 1j*eta*Green_Gamma_width)
        D_n2k_4 = complex(1.,0.)/(xi - xi_n2k + hw + 1j*eta*Green_Gamma_width)

        f_xi_ehw = function_fermi_occupation(list_xi-eta*hw,0.,beta)[0]

        Term_1  =  ( f_xi - f_xi_ehw ) * (D_n2k_1 - D_n2k_2) 

        Term_2  =   f_xi * (D_n2k_3 - D_n2k_4) 

        common_factor = prefactor*kernel*D_n1k*D_n3kq 

        J_integrand_df += common_factor*Term_1  
        J_integrand_f  += common_factor*Term_2  

    J_integrand = J_integrand_df+ J_integrand_f  

    return J_integrand

def get_approximate_J_df_integrand(xi, xi_n1k, xi_n2k, xi_n3kq, hw, mu, beta, kernel_Gamma_width, Green_Gamma_width):
    """
    This function computes the exact integrand leading to the J function.
    """

    J_integrand_df = complex(0.,0.)
    J_integrand_f  = complex(0.,0.)

    list_xi = N.array([xi])

    KR = get_KR(list_xi+mu,kernel_Gamma_width)[0]
    KI = get_KI(list_xi+mu,kernel_Gamma_width)[0]

    df = d_Fermi_dxi(list_xi,0.,beta)[0]

    for eta in [-1,1]:
        prefactor = eta/(2.*N.pi*1j)

        kernel = KR-1j*eta*KI

        D_n1k  = complex(1.,0.)/(xi - xi_n1k  + 1j*eta*Green_Gamma_width)
        D_n3kq = complex(1.,0.)/(xi - xi_n3kq + 1j*eta*Green_Gamma_width)

        D_n2k_1 = complex(1.,0.)/(xi - xi_n2k -eta*hw - 1j*Green_Gamma_width)
        D_n2k_2 = complex(1.,0.)/(xi - xi_n2k -eta*hw + 1j*Green_Gamma_width)

        D_n2k_3 = complex(1.,0.)/(xi - xi_n2k - hw + 1j*eta*Green_Gamma_width)
        D_n2k_4 = complex(1.,0.)/(xi - xi_n2k + hw + 1j*eta*Green_Gamma_width)


        Term_1  =  eta*hw*df* (D_n2k_1 - D_n2k_2) 

        common_factor = prefactor*kernel*D_n1k*D_n3kq 

        J_integrand_df += common_factor*Term_1  


    return J_integrand_df

def get_J_df_approx(xi_n1k, xi_n2k, xi_n3kq, hw, mu, beta, kernel_Gamma_width, Green_Gamma_width):
    """
    This function computes the approximate form for the df term of J, which is no longer an integrand.
    """

    J_integrand_df = complex(0.,0.)
    J_integrand_f  = complex(0.,0.)

    xi = 0.

    list_xi = N.array([xi])

    KR = get_KR(list_xi+mu,kernel_Gamma_width)[0]
    KI = get_KI(list_xi+mu,kernel_Gamma_width)[0]


    for eta in [-1,1]:
        prefactor = eta/(2.*N.pi*1j)

        kernel = KR-1j*eta*KI

        D_n1k  = complex(1.,0.)/(xi - xi_n1k  + 1j*eta*Green_Gamma_width)
        D_n3kq = complex(1.,0.)/(xi - xi_n3kq + 1j*eta*Green_Gamma_width)

        D_n2k_1 = complex(1.,0.)/(xi - xi_n2k -eta*hw - 1j*Green_Gamma_width)
        D_n2k_2 = complex(1.,0.)/(xi - xi_n2k -eta*hw + 1j*Green_Gamma_width)

        D_n2k_3 = complex(1.,0.)/(xi - xi_n2k - hw + 1j*eta*Green_Gamma_width)
        D_n2k_4 = complex(1.,0.)/(xi - xi_n2k + hw + 1j*eta*Green_Gamma_width)


        Term_1  =  ( -eta * hw ) * (D_n2k_1 - D_n2k_2) 


        common_factor = prefactor*kernel*D_n1k*D_n3kq 

        J_integrand_df += common_factor*Term_1  


    return J_integrand_df



