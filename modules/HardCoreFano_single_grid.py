#!/usr/local/bin/python

#================================================================================
#
# This program gets parameters from input and computes the HCF spectrum.
#
#================================================================================

# import  modules
import sys
import os
import socket


hostname = socket.gethostname().split('.')[0]

if 'mbpdeinllation' in hostname:
    top_dir = '/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_10.0/'
else:
    top_dir = '/RQusagers/roussea4/python_depository/HCF_10.0/'


sys.path.append(top_dir)

import module_Driver
reload(module_Driver)
from module_Driver import *

#----------------------------------------
# Read input parameters from command line
#----------------------------------------
input_error   = False

args = sys.argv[1:]
if len(args) != 21:
        input_error = True

try:
        mu     = N.float(args[0])
        T      = N.float(args[1])

        nmax_fine   = N.int(args[2])


        Green_Gamma_width   = N.float(args[3])
        kernel_Gamma_width  = N.float(args[4])


        hw_ph  = N.float(args[5])
        q_ph1  = N.float(args[6])
        q_ph2  = N.float(args[7])

        re_E1x = N.float(args[8])
        im_E1x = N.float(args[9])
        re_E1y = N.float(args[10])
        im_E1y = N.float(args[11])
        re_E1z = N.float(args[12])
        im_E1z = N.float(args[13])

        re_E2x = N.float(args[14])
        im_E2x = N.float(args[15])
        re_E2y = N.float(args[16])
        im_E2y = N.float(args[17])
        re_E2z = N.float(args[18])
        im_E2z = N.float(args[19])
        filename = args[20]
except:
        input_error   = True


if input_error:
        print 'Something is wrong with input parameters!'
        print 'HardCoreFano.py  <mu> <T>  <nk_grid_fine>  <Green_Gamma_width> <kernel_Gamma_width> <hw_ph> <q_ph> <E_ph> <output_filename>'
        print '                mu            :     chemical potential, in eV'
        print '                T             :     temperature, in Kelvin'
        print '                nk_grid_fine  :     parameter specifying how dense the k-grid will be '
        print '           Green_Gamma_width  :     width for Green functions, in eV'
        print '           kernel_Gamma_width :     Broadening used in the scattering kernel function, in eV'
        print '                hw_ph         :     phonon energy, in eV'
        print '                q_ph          :     phonon q vector; two real numbers'
        print '                E_ph          :     phonon polarization vector; 12 real numbers (alternating real and imaginary part)'
        print '        output_filename       :     name of netcdf file where to write the data'
        sys.exit()


# Setup parameters

beta = 1./(kB*T)

q_ph = N.array([q_ph1  ,q_ph2])  

E_ph = N.array([ re_E1x +1j*im_E1x , re_E1y +1j*im_E1y , re_E1z +1j*im_E1z ,
                 re_E2x +1j*im_E2x , re_E2y +1j*im_E2y , re_E2z +1j*im_E2z ])


#================================================================================
# Generate data!
#================================================================================

grid = TesselationGrid(nmax_fine)

#================================================================================
#       Compute the Hq function 
#================================================================================

OkComputer = Compute_Loop_Function(mu, beta, q_ph, E_ph, hw_ph, grid, kernel_Gamma_width, Green_Gamma_width)

OkComputer.Compute_R_and_I()
nmax_coarse = 0
nblock = 0
write_to_file(OkComputer, nmax_coarse, nmax_fine, nblock, hw_ph, filename)
