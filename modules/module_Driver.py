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

import module_I
reload(module_I)
from module_I              import *

import module_Tests
reload(module_Tests)
from module_Tests          import *


import module_Integrators
reload(module_Integrators)
from module_Integrators    import *


import module_Compute_Loop_Function
reload(module_Compute_Loop_Function)
from module_Compute_Loop_Function import *

import module_Kramers_Kronig
reload(module_Kramers_Kronig)
from module_Kramers_Kronig import *

import module_NETCDF
reload(module_NETCDF)
from module_NETCDF         import *


#================================================================================
# useful commands
#
#================================================================================
def build_command(executable,mu,T,nk_max_coarse, nk_max_fine, nblocks,n_hw,hw_max,Gamma,hw_ph,q_ph,E_ph,output_filename):

    command = [executable,
               '%8.4f'%mu,
               '%8.2f'%T,
                  '%i'%nk_max_coarse,
                  '%i'%nk_max_fine,
                  '%i'%nblocks,
                  '%i'%n_hw,
               '%8.4f'%hw_max,
               '%8.4f'%Gamma,
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


def build_command_imaginary(executable,mu,T,n_hw,hw_max,Gamma,hw_ph,q_ph,E_ph,output_filename):

    command = [executable,
               '%8.4f'%mu,
               '%8.2f'%T,
                  '%i'%n_hw,
               '%8.4f'%hw_max,
               '%8.4f'%Gamma,
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

def build_command_matsubara_sum(executable,mu,T,nk_max_coarse, nk_max_fine, nblocks,n_hw,hw_max,Gamma, \
                        matsubara_cutoff,hw_ph,q_ph,E_ph,output_filename):

    command = [executable,
               '%8.4f'%mu,
               '%8.2f'%T,
                  '%i'%nk_max_coarse,
                  '%i'%nk_max_fine,
                  '%i'%nblocks,
                  '%i'%n_hw,
               '%8.4f'%hw_max,
               '%8.4f'%Gamma,
               '%8.4f'%matsubara_cutoff,
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

def build_string(mu,T,nk_max_coarse, nk_max_fine, nblocks,n_hw,hw_max,Gamma, hw_ph,q_ph,E_ph,output_filename):

    string =   '%8.4f '%mu+\
               '%8.2f '%T+\
                  '%i '%nk_max_coarse+\
                  '%i '%nk_max_fine+\
                  '%i '%nblocks+\
                  '%i '%n_hw+\
               '%8.4f '%hw_max+\
               '%8.4f '%Gamma+\
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


def build_string_matsubara_sum(mu,T,nk_max_coarse, nk_max_fine, nblocks,n_hw,hw_max,Gamma,matsubara_cutoff, hw_ph,q_ph,E_ph,output_filename):

    string =   '%8.4f '%mu+\
               '%8.2f '%T+\
                  '%i '%nk_max_coarse+\
                  '%i '%nk_max_fine+\
                  '%i '%nblocks+\
                  '%i '%n_hw+\
               '%8.4f '%hw_max+\
               '%8.4f '%Gamma+\
               '%8.4f '%matsubara_cutoff+\
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


