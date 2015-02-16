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
    top_dir = '/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_5.0/'
else:
    top_dir = '/RQusagers/roussea4/python_depository/HCF_5.0/'


sys.path.append(top_dir)

import module_Driver
reload(module_Driver)
from module_Driver import *

#----------------------------------------
# Read input parameters from command line
#----------------------------------------
input_error   = False

args = sys.argv[1:]
if len(args) != 24:
        input_error = True

try:
        mu     = N.float(args[0])
        T      = N.float(args[1])

        nmax_coarse  = N.int(args[2])
        nmax_fine    = N.int(args[3])
        nblock       = N.int(args[4])

        n_hw                = N.int(args[5])
        hw_max              = N.float(args[6])
        delta_width         = N.float(args[7])

        hw_ph  = N.float(args[8])
        q_ph1  = N.float(args[9])
        q_ph2  = N.float(args[10])

        re_E1x = N.float(args[11])
        im_E1x = N.float(args[12])
        re_E1y = N.float(args[13])
        im_E1y = N.float(args[14])
        re_E1z = N.float(args[15])
        im_E1z = N.float(args[16])

        re_E2x = N.float(args[17])
        im_E2x = N.float(args[18])
        re_E2y = N.float(args[19])
        im_E2y = N.float(args[20])
        re_E2z = N.float(args[21])
        im_E2z = N.float(args[22])
        filename = args[23]
except:
        input_error   = True


if input_error:
        print 'Something is wrong with input parameters!'
        print 'HardCoreFano_weak_scattering.py  <mu> <T> <nk_grid_coarse> <nk_grid_fine> <nblocks> <n_hw> <hw_max> <hw_ph> <q_ph> <E_ph> <output_filename>'
        print '                mu            :     chemical potential, in eV'
        print '                T             :     temperature, in Kelvin'
        print '                nk_grid_coarse:     parameter specifying how dense the coarse k-grid will be for the smooth grid'
        print '                nk_grid_fine  :     parameter specifying how dense the fine k-grid will be for the smooth grid'
        print '                nblocks       :     parameter specifying how many coarse blocks will be made fine for the smooth grid'
        print '                n_hw          :     number of points on frequency grid'
        print '                hw_max        :     maximum frequency on frequency grid, in eV'
        print '                delta_width   :     width for delta-functions, in eV'
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
include_Gamma = False

grid          = TesselationDoubleGrid(nmax_coarse, nmax_fine, nblock, include_Gamma,clip_grid=False,clip_energy=0.)

# We suppose that the grid is regular and goes from 100 meV to hw_max.

hw_min = 0.100 # meV

d_hw         = (hw_max-hw_min)/(n_hw-1.)
iloop        = N.arange(n_hw)
list_hw      = hw_min+d_hw*iloop


#================================================================================
#       Compute the Hq function 
#================================================================================

OkComputer = Compute_Loop_Function_Product(mu, beta, q_ph, E_ph, hw_ph, grid, list_hw, delta_width)

OkComputer.Compute_Hq_Product()

write_to_file_weak_coupling(OkComputer, nmax_coarse, nmax_fine, nblock, hw_ph, filename) 
