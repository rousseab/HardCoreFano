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
if len(args) != 28:
        input_error = True

try:
        mu     = N.float(args[0])
        T      = N.float(args[1])

        nmax_coarse_smooth   = N.int(args[2])
        nmax_fine_smooth     = N.int(args[3])
        nblock_smooth        = N.int(args[4])

        nmax_coarse_singular = N.int(args[5])
        nmax_fine_singular   = N.int(args[6])
        nblock_singular      = N.int(args[7])



        n_hw                = N.int(args[8])
        hw_max              = N.float(args[9])
        delta_width         = N.float(args[10])
        kernel_Gamma_width  = N.float(args[11])


        hw_ph  = N.float(args[12])
        q_ph1  = N.float(args[13])
        q_ph2  = N.float(args[14])

        re_E1x = N.float(args[15])
        im_E1x = N.float(args[16])
        re_E1y = N.float(args[17])
        im_E1y = N.float(args[18])
        re_E1z = N.float(args[19])
        im_E1z = N.float(args[20])

        re_E2x = N.float(args[21])
        im_E2x = N.float(args[22])
        re_E2y = N.float(args[23])
        im_E2y = N.float(args[24])
        re_E2z = N.float(args[25])
        im_E2z = N.float(args[26])
        filename = args[27]

except:
        input_error   = True


if input_error:
        print 'Something is wrong with input parameters!'
        print 'HardCoreFano.py  <mu> <T> <nk_grid_coarse_smooth> <nk_grid_fine_smooth> <nblocks_smooth> <<nk_grid_coarse_singular> <nk_grid_fine_singular> <nblocks_singular> <n_hw> <hw_max> <Gamma> <hw_ph> <q_ph> <E_ph> <output_filename>'
        print '                mu            :     chemical potential, in eV'
        print '                T             :     temperature, in Kelvin'
        print ' SMOOTH         nk_grid_coarse:     parameter specifying how dense the coarse k-grid will be for the smooth grid'
        print ' SMOOTH         nk_grid_fine  :     parameter specifying how dense the fine k-grid will be for the smooth grid'
        print ' SMOOTH         nblocks       :     parameter specifying how many coarse blocks will be made fine for the smooth grid'
        print ' SINGULAR       nk_grid_coarse:     parameter specifying how dense the coarse k-grid will be for the singular grid'
        print ' SINGULAR       nk_grid_fine  :     parameter specifying how dense the fine k-grid will be for the singular grid'
        print ' SINGULAR       nblocks       :     parameter specifying how many coarse blocks will be made fine for the singular grid'
        print '                n_hw          :     number of points on frequency grid'
        print '                hw_max        :     maximum frequency on frequency grid, in eV'
        print '                delta_width   :     width for delta-functions, in eV'
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
include_Gamma = False

grid_smooth = TesselationDoubleGrid(nmax_coarse_smooth, nmax_fine_smooth, nblock_smooth, include_Gamma,clip_grid=False,clip_energy=0.)

clip_energy   = N.abs(mu)+0.5
grid_singular = TesselationDoubleGrid(nmax_coarse_singular, nmax_fine_singular, nblock_singular, include_Gamma,clip_grid=True,clip_energy=clip_energy)

# We suppose that the grid is regular and goes from 100 meV to hw_max.

hw_min = 0.100 # meV

d_hw         = (hw_max-hw_min)/(n_hw-1.)
iloop        = N.arange(n_hw)
list_hw      = hw_min+d_hw*iloop


#================================================================================
#       Compute the Hq function 
#================================================================================

OkComputer = Compute_Loop_Function_Product(type, mu, beta, q_ph, E_ph, hw_ph, grid_singular, grid_smooth, list_hw, kernel_Gamma_width, delta_width)

OkComputer.Compute_Hq_Product()

write_to_file(OkComputer, nmax_coarse, nmax_fine, nblock, hw_ph, filename)
