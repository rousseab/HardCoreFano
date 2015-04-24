#================================================================================
#
# This test computes the integral for J numerically, and plots the 
# results for various values of the parameters xi_1, xi_3, hw.
#
#================================================================================
import common
reload(common)
from common import *

import module_numerical_tests 
reload(module_numerical_tests )
from module_numerical_tests import *


import numpy as N
import matplotlib.pyplot as plt


#parameters 
T = 300.
beta = 1./(kB*T)

kernel_Gamma_width = 0.200 # eV
Green_Gamma_width  = 0.025
mu = -0.400
hw = 0.1755


Integrand_generator = NumericalJ( hw, mu, beta, kernel_Gamma_width, Green_Gamma_width)

global_integrand = Integrand_generator.get_exact_J_integrand

list_xi_n3kq = N.arange(-1.0,1.01,0.025)
list_xi_n1k  = N.arange(-1.0,1.01,0.025)



filename = 'J_Green_Gamma=%i_meV_mu=%i_meV_hw=%4.1f_meV.nc'%(1000*Green_Gamma_width,1000*mu,1000*hw)
handler = NetcdfHandler(filename)
handler.open_ncfile(mode='w')
handler.write_attributes(mu, beta, Green_Gamma_width, kernel_Gamma_width, hw, list_xi_n1k, list_xi_n3kq)

for state_n2 in ['intraband','interband']:

    if state_n2 == 'intraband':
        list_xi_n2k = list_xi_n1k 
    elif state_n2 == 'interband':
        list_xi_n2k = -list_xi_n1k -2*mu

    J_values = N.complex(0.,0.)*N.zeros([list_xi_n1k.shape[0],list_xi_n3kq.shape[0]])

    for ik, (xi_n1k,xi_n2k) in enumerate(zip(list_xi_n1k, list_xi_n2k)):

        global_parameters = {'xi_n1k':xi_n1k,'xi_n2k':xi_n2k}

        list_re_J_numeric, list_im_J_numeric = \
                Drive_Numerical_Integrals(list_xi_n3kq, global_integrand, 'xi_n3kq', global_parameters)

        J_values[ik,:] = complex(1.,0)*list_re_J_numeric+ complex(0.,1.)*list_im_J_numeric


    handler.write_J_array(J_values,'numerical_J_'+state_n2)

handler.close_ncfile()


