#================================================================================
#
# Generate some data
#
#
#
#================================================================================

import time
import warnings
warnings.filterwarnings("ignore")

import sys
import os
import numpy as N

top_dir = '/RQusagers/roussea4/python_depository/HCF_5.0/'

sys.path.append(top_dir+'modules/')

from module_Driver import *

submit_file_tmpl ="""#!/bin/bash
#PBS -N grph_%i
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=1GB

module load intel-compilers/12.0.4.191
module load python/2.7.2
module load numpy/1.7.1
module load netcdf/4.1.3
module load scipy/0.12.0python2.7


export MV2_ENABLE_AFFINITY=0 #To avoid node problems
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

ulimit -s unlimited   # unlimited stack size!


PYTHONPATH=/RQusagers/roussea4/python_depository/lib/python2.7/site-packages:$PYTHONPATH

python=/home/apps/Logiciels/Python/python-2.7.2/bin/python

cd ${PBS_O_WORKDIR}   # go to work directory

HCK="$python /RQusagers/roussea4/python_depository/HCF_5.0/modules/HardCoreFano.py"
"""
#----------------------------------------
# Global parameters
#----------------------------------------
T    = 300.  # Kelvin
kB = 8.6173324e-5 # eV/K
beta = 1./(kB*T)                                                                                                                                                            


SUBMIT = True
#SUBMIT = False

list_nu_index = [5,4,3,2]

mu = -0.400

kernel_Gamma_width = 0.050
delta_width        = 0.050


n_hw   = 251
hw_max = 0.250

k_nmax_coarse_singular = 8
k_nmax_fine_singular   = 256
k_n_blocks_coarse_to_fine_singular = 2
kg_singular = [k_nmax_coarse_singular, k_nmax_fine_singular, k_n_blocks_coarse_to_fine_singular]

k_nmax_coarse_smooth = 32
k_nmax_fine_smooth   = 128
k_n_blocks_coarse_to_fine_smooth = 8
kg_smooth = [k_nmax_coarse_smooth, k_nmax_fine_smooth, k_n_blocks_coarse_to_fine_smooth]


q_nmax_coarse = 32
q_nmax_fine   = 128
q_n_blocks_coarse_to_fine = 5

qg = [q_nmax_coarse ,q_nmax_fine , q_n_blocks_coarse_to_fine]

q_include_Gamma = True

submit_max    = 70

#================================================================================
# Phonon q point grid
#================================================================================

mauri_filepath   = top_dir+'modules/' 


list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()
Renormalizor     = CompteMauriRenormalize(mauri_filepath,list_FC,list_R)

tol    = 1e-10

Q_grid = TesselationDoubleGrid(q_nmax_coarse, q_nmax_fine, q_n_blocks_coarse_to_fine,q_include_Gamma )


Q_wedge = Q_grid.list_wedges[0]
list_q  = []


nk = len(Q_wedge.list_k)
for i,k in enumerate(Q_wedge.list_k):
    print 'checking vector %i of %i'%(i+1,nk)

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



number_of_q_vectors = len(list_q )

#================================================================================
# Generate data!
#================================================================================


filename_template  = 'HCF_iq=%i_nu=%i.nc'

counter  = 0



workdir = 'delta=%i_meV_Gamma=%i_meV_mu=%i_meV'%(1000.*delta_width,1000.*kernel_Gamma_width,1000.*mu)+\
          '_nq=%i_%i_%i'%(qg[0],qg[1],qg[2])+\
          '_nk_smooth=%i_%i_%i'%(kg_smooth[0],kg_smooth[1],kg_smooth[2])+\
          '_nk_singular=%i_%i_%i/'%(kg_singular[0],kg_singular[1],kg_singular[2])

print 'generating %s...'%workdir 

#os.mkdir(workdir)
os.chdir(workdir)


# make the scattering kernel
"""
SK = ScatteringKernel(mu,beta,kernel_Gamma_width)
kernel_filename ='scattering_spline.nc'
SK.build_scattering_kernel()
SK.build_and_write_spline_parameters(kernel_filename)
"""


# first, generate all the independent tasks to be performed.

list_task_strings = []
iq_index = 0
for q_ph in  list_q:

    iq_index += 1
    list_hw_ph_NOT_RENORMALIZED, list_hw_ph, list_Eph = Renormalizor.get_frequency_and_polarization(q_ph)

    for nu_index in list_nu_index:
        E_ph  = list_Eph[:,nu_index]
        hw_ph = list_hw_ph[nu_index]


        filename  = filename_template%(iq_index,nu_index)

        string = build_string(mu,T,kg_smooth[0],kg_smooth[1],kg_smooth[2],\
                                kg_singular[0],kg_singular[1],kg_singular[2],\
                                n_hw,hw_max,delta_width, kernel_Gamma_width, hw_ph,q_ph,E_ph,filename)

        list_task_strings.append('$HCK  %s'%string)


# Find out how many tasks should go in every submit script, to avoid having more than submit_max scripts
jobs_per_script = N.ceil(1.0*len(list_task_strings)/submit_max)


# fill the submit_scripts

for i, task in enumerate(list_task_strings):

    if N.mod(i,jobs_per_script ) == 0:
        try:
            file.close()
            if SUBMIT == True:
                os.system('qsub %s'%submit_filename)
        except:
            pass

        counter += 1
        submit_filename = 'submit_%i.pbs'%counter
        file = open(submit_filename,'w')
        print >> file, submit_file_tmpl%(counter)


    print >> file, task


file.close()
if SUBMIT == True:
    os.system('qsub %s'%submit_filename)


os.chdir('../')

