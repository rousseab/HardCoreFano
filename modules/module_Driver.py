#================================================================================
#
#           Module Driver
#           ==========================================
#
#       This module imports all other modules and makes running jobs possible.
#
#
#================================================================================

#----------------------------------------
# import modules
#----------------------------------------

import module_Constants       
reload(module_Constants)
from module_Constants      import *

import module_Path
reload(module_Path)
from module_Path           import *

import module_Functions
reload(module_Functions)
from module_Functions      import *

import  module_MatrixList
reload(module_MatrixList)
from module_MatrixList import *


import module_Phonons
reload(module_Phonons)
from module_Phonons        import *

import module_ForceConstants_Dubay
reload(module_ForceConstants_Dubay)
from module_ForceConstants_Dubay import *


import module_Grids           
reload(module_Grids)
from module_Grids          import *

import module_D6h
reload(module_D6h)
from module_D6h            import *

import module_M
reload(module_M)
from module_M              import *

import module_HardCoreKernel  
reload(module_HardCoreKernel)
from module_HardCoreKernel import *

import module_J
reload(module_J)
from module_J              import *

import module_Tests
reload(module_Tests)
from module_Tests          import *


import module_Integrators
reload(module_Integrators)
from module_Integrators    import *


import module_Kramers_Kronig
reload(module_Kramers_Kronig)
from module_Kramers_Kronig import *

import module_Compute_Loop_Function
reload(module_Compute_Loop_Function)
from module_Compute_Loop_Function import *



import module_NETCDF
reload(module_NETCDF)
from module_NETCDF         import *


#================================================================================
# useful commands
#
#================================================================================
def build_command(executable,mu,T,nk_max_coarse, nk_max_fine, nblocks,\
                        Green_Gamma_width, kernel_Gamma_width,hw_ph,q_ph,E_ph,output_filename):

    command = [executable,
               '%8.4f'%mu,
               '%8.2f'%T,
                  '%i'%nk_max_coarse,
                  '%i'%nk_max_fine,
                  '%i'%nblocks,
               '%8.4f'%Green_Gamma_width,
               '%8.4f'%kernel_Gamma_width,
             '%16.12f'%hw_ph,
             '%20.16f'%q_ph[0],
             '%20.16f'%q_ph[1],
             '%20.16f'%N.real(E_ph[0]),
             '%20.16f'%N.imag(E_ph[0]),
             '%20.16f'%N.real(E_ph[1]),
             '%20.16f'%N.imag(E_ph[1]),
             '%20.16f'%N.real(E_ph[2]),
             '%20.16f'%N.imag(E_ph[2]),
             '%20.16f'%N.real(E_ph[3]),
             '%20.16f'%N.imag(E_ph[3]),
             '%20.16f'%N.real(E_ph[4]),
             '%20.16f'%N.imag(E_ph[4]),
             '%20.16f'%N.real(E_ph[5]),
             '%20.16f'%N.imag(E_ph[5]),
                  '%s'%output_filename]

    return command

def build_string(mu,T,nk_max_coarse, nk_max_fine, nblocks,\
                    Green_Gamma_width,kernel_Gamma_width,hw_ph,q_ph,E_ph,output_filename):

    string =   '%8.4f '%mu+\
               '%8.2f '%T+\
                  '%i '%nk_max_coarse+\
                  '%i '%nk_max_fine+\
                  '%i '%nblocks+\
               '%8.4f '%Green_Gamma_width+\
               '%8.4f '%kernel_Gamma_width+\
             '%16.12f '%hw_ph+\
             '%20.16f '%q_ph[0]+\
             '%20.16f '%q_ph[1]+\
             '%20.16f '%N.real(E_ph[0])+\
             '%20.16f '%N.imag(E_ph[0])+\
             '%20.16f '%N.real(E_ph[1])+\
             '%20.16f '%N.imag(E_ph[1])+\
             '%20.16f '%N.real(E_ph[2])+\
             '%20.16f '%N.imag(E_ph[2])+\
             '%20.16f '%N.real(E_ph[3])+\
             '%20.16f '%N.imag(E_ph[3])+\
             '%20.16f '%N.real(E_ph[4])+\
             '%20.16f '%N.imag(E_ph[4])+\
             '%20.16f '%N.real(E_ph[5])+\
             '%20.16f '%N.imag(E_ph[5])+\
                  '%s '%output_filename

    return string


def Drive_Numerical_Integrals(list_x_numeric, global_integrand, x_str_name, global_parameters):
    """
    This routine performs numerical integrals for a range of parameters.

    We assume that the integrand is of the form
        global_integrand(xi, x, --- global_parameters ---)
    """
    from functools import  partial
    from scipy.integrate import quad

    print ' ==== Performing Numerical Integration ====='

    list_re_I_numeric = []
    list_im_I_numeric = []

    parameters = deepcopy(global_parameters)
    for x in list_x_numeric:

        # update the parameters dictionary
        parameters[x_str_name] = x
        integrand = partial(global_integrand, **parameters)

        re_f = lambda x: N.real(integrand(x))
        im_f = lambda x: N.imag(integrand(x))

        re_I, rerr = quad(re_f, -N.infty, N.infty)
        im_I, ierr = quad(im_f, -N.infty, N.infty)
        
        print '   doing x = %4.3f: '%x 
        print '             - Rerr = %8.4e'%rerr
        print '             - Ierr = %8.4e'%ierr

        list_re_I_numeric.append(re_I)
        list_im_I_numeric.append(im_I)

    list_re_I_numeric = N.array(list_re_I_numeric )
    list_im_I_numeric = N.array(list_im_I_numeric )

    return list_re_I_numeric, list_im_I_numeric 

