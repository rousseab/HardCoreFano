#================================================================================
#
# Generate some data
#
#================================================================================
import sys
import os
import time
import subprocess as SUB

top_dir = '/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_5.0/'

sys.path.append(top_dir+'modules/')

import module_Driver
reload(module_Driver)
from module_Driver import *

#executable
path= top_dir+'modules/'
HCK = path+'HardCoreFano.py'


#----------------------------------------
# Global parameters
#----------------------------------------

max_processes = 8

mu   =-0.400 # eV
T    = 300.  # Kelvin

kB   = 8.6173324e-5 # eV/K
beta = 1./(kB*T)

list_nu_index = [2,3,4,5]

nkmax_coarse_singular =  8
nkmax_fine_singular   = 32
nk_blocks_coarse_to_fine_singular = 3

nkmax_coarse_smooth =  8
nkmax_fine_smooth = 16
nk_blocks_coarse_to_fine_smooth = 2


nqmax_coarse =  8
nqmax_fine   = 32
nq_blocks_coarse_to_fine = 3

kernel_Gamma_width  = 0.200 # meV
delta_width         = 0.050 # meV

n_hw   = 151
hw_max = 0.250


#================================================================================
# Phonon q point grid
#================================================================================


list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()

mauri_filepath= path

Renormalizor = CompteMauriRenormalize(mauri_filepath,list_FC,list_R)

q_include_Gamma = True

Q_grid = TesselationDoubleGrid(nqmax_coarse, nqmax_fine, nq_blocks_coarse_to_fine,q_include_Gamma )



Q_wedge = Q_grid.list_wedges[0]
list_q  = []

tol    = 1e-10

for k in Q_wedge.list_k:

    if len(list_q) == 0:
        new = True
    else:
        list_dq = N.array(list_q)-k

        distances = N.sqrt(N.sum(list_dq**2 ,axis=1))

        if distances.min() < tol:
            new = False
        else:
            new = True

    if new:
        list_q.append(k)

list_q  = N.array(list_q)

number_of_q_vectors = len(list_q )


x = Q_wedge.list_k[:,0]
y = Q_wedge.list_k[:,1]

dx = x[:,N.newaxis]-x[N.newaxis,:]
dy = y[:,N.newaxis]-y[N.newaxis,:]

dz = N.sqrt(dx**2+dy**2)

tol = 1e-10
I,J  = N.where(dz < tol)

repeated_indices = N.where(I!=J)[0]

i = I[repeated_indices] 
j = J[repeated_indices] 

k1 = Q_wedge.list_k[i]
k2 = Q_wedge.list_k[j]

