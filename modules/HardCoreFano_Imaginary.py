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

if hostname == 'ferron':
        top_dir = '/Users/Shared/Bruno/work/projects/Graphene_Fano/python_work/HardCoreKernelFano/'
elif hostname == 'briaree1':
        top_dir = '/RQusagers/roussea4/python_depository/HCF/'
else:
    top_dir = '/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_4.0/modules/'


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

        n_hw   = N.int(args[2])
        hw_max = N.float(args[3])
        Gamma  = N.float(args[4])


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
        print 'HardCoreFano_Imaginary.py <mu> <T> <n_hw> <hw_max> <Gamma> <hw_ph> <q_ph> <E_ph> <output_filename>'
        print '                mu            :     chemical potential, in eV'
        print '                T             :     temperature, in Kelvin'
        print '                n_hw          :     number of points on frequency grid'
        print '                hw_max        :     maximum frequency on frequency grid, in eV'
        print '                Gamma         :     width, in eV, for the scattering kernel'
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
# We suppose that the grid is regular and goes from zero to hw_max.

d_hw         = hw_max/(n_hw-1.)
iloop        = N.arange(n_hw)
list_hw      = d_hw*iloop


#================================================================================
#       Compute the Hq function 
#================================================================================

OkComputer = Compute_Imaginary_Loop_Function( mu, beta, q_ph, E_ph, hw_ph,  list_hw, Gamma)

OkComputer.Compute_Hq()

write_to_file(OkComputer, 0, 0, 0, hw_ph, filename)
