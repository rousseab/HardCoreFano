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

max_processes = 4

mu   =-0.400 # eV
T    = 300.  # Kelvin

kB   = 8.6173324e-5 # eV/K
beta = 1./(kB*T)

list_nu_index = [2,3,4,5]

nkmax_coarse =  8
nkmax_fine   = 32
nk_blocks_coarse_to_fine = 3

nqmax_coarse =  8
nqmax_fine   = 32
nq_blocks_coarse_to_fine = 3

kernel_Gamma_width  = 0.050 # meV
delta_width         = 0.050 # meV

n_hw   = 251
hw_max = 0.250

type_of_integral = 'smooth'

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

    new = True

    list_dq = N.array(list_q)-k

    distances = N.sqrt(N.sum(list_dq**2 ,axis=1))

    if distances.min() < tol:
        new = False

    if new:
        list_q.append(k)

list_q  = N.array(list_q)

number_of_q_vectors = len(list_q )

#================================================================================
# Generate data!
#================================================================================

filename_template  = 'HCF_nq=%i_%i_%i_nk=%i_%i_%i_iq=%i_nu=%i.nc'

log = open('calculation.log','w')


set_of_processes = set()

work_dir = 'nq=%i_%i_%i_nk=%i_%i_%i_mu=%4.3f/'%(nqmax_coarse,nqmax_fine,nq_blocks_coarse_to_fine, nkmax_coarse, nkmax_fine,nk_blocks_coarse_to_fine,mu)
os.mkdir(work_dir )
os.chdir(work_dir )

# make the scattering kernel
SK = ScatteringKernel(mu,beta,kernel_Gamma_width)
kernel_filename ='scattering_spline.nc'
SK.build_scattering_kernel()
SK.build_and_write_spline_parameters(kernel_filename)


for nu_index in list_nu_index:
    iq_index = 0
    for q_ph in  list_q:

        iq_index  += 1
        list_hw_ph_NOT_RENORMALIZED, list_hw_ph, list_Eph= Renormalizor.get_frequency_and_polarization(q_ph)

        hw_ph = list_hw_ph[nu_index]
        E_ph  = list_Eph[:,nu_index]

        filename = filename_template%(nqmax_coarse,nqmax_fine,nq_blocks_coarse_to_fine, nkmax_coarse, nkmax_fine,\
                                        nk_blocks_coarse_to_fine,iq_index,nu_index)
        
        command  = build_command(HCK,type_of_integral,mu,T,nkmax_coarse,nkmax_fine, nk_blocks_coarse_to_fine, n_hw,hw_max,\
                                    delta_width,kernel_Gamma_width,hw_ph,q_ph1,E_ph1,filename_plus)

        # take out finished jobs
        while len(set_of_processes) >= max_processes:

            time.sleep(2)
            set_of_finished_processes = set()
            for job in set_of_processes:
                if job.poll() is not None:
                    set_of_finished_processes.add(job)

            set_of_processes.difference_update(set_of_finished_processes)


        print 'Doing nu =%i, iq = %i of %i q-points'%(nu_index,iq_index,number_of_q_vectors)
        print >> log,'Doing nu =%i, iq = %i of %i q-points'%(nu_index,iq_index,number_of_q_vectors)
        log.flush()

        job = SUB.Popen(command)
        set_of_processes.add(job)


os.chdir('../')

log.close()
