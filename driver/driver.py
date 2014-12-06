#================================================================================
#
# Generate some data
#
#================================================================================
import sys
import os
import time
import subprocess as SUB

top_dir = '/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_4.0/'

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

Gamma         = 0.050 # meV

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


tol    = 1e-10


list_plus_q  = []
list_minus_q = []
for wedge in  Q_grid.list_wedges:

        for k in wedge.list_k:

                new = True

                for q1, q2 in zip(list_plus_q, list_minus_q):

                        test1 = N.linalg.norm(k-q1) < tol
                        test2 = N.linalg.norm(k-q2) < tol
                        if test1 or test2:
                                new = False
                                break

                if new:
                        list_plus_q.append(k)
                        list_minus_q.append(-k)

list_plus_q  = N.array(list_plus_q)
list_minus_q = N.array(list_minus_q)


number_of_q_pairs = len(list_plus_q )

#================================================================================
# Generate data!
#================================================================================

filename_plus_template  = 'HCF_nq=%i_%i_%i_nk=%i_%i_%i_iq_plus=%i_nu=%i.nc'
filename_minus_template = 'HCF_nq=%i_%i_%i_nk=%i_%i_%i_iq_minus=%i_nu=%i.nc'


log = open('calculation.log','w')


set_of_processes = set()

work_dir = 'nq=%i_%i_%i_nk=%i_%i_%i_mu=%4.3f/'%(nqmax_coarse,nqmax_fine,nq_blocks_coarse_to_fine, nkmax_coarse, nkmax_fine,nk_blocks_coarse_to_fine,mu)
os.mkdir(work_dir )
os.chdir(work_dir )

# make the scattering kernel
SK = ScatteringKernel(mu,beta,Gamma)
kernel_filename ='scattering_spline.nc'
SK.build_scattering_kernel()
SK.build_and_write_spline_parameters(kernel_filename)


for nu_index in list_nu_index:
    iq_index = 0
    for q_ph1, q_ph2 in  zip(list_plus_q,  list_minus_q):

        iq_index  += 1

        list_hw_ph_NOT_RENORMALIZED, list_hw_ph, list_Eph= Renormalizor.get_frequency_and_polarization(q_ph1)


        hw_ph = list_hw_ph[nu_index]

        E_ph1 = list_Eph[:,nu_index]
        E_ph2 = N.conjugate(E_ph1)


        filename_plus  = filename_plus_template%(nqmax_coarse,nqmax_fine,nq_blocks_coarse_to_fine, nkmax_coarse, nkmax_fine,nk_blocks_coarse_to_fine,iq_index,nu_index)
        filename_minus = filename_minus_template%(nqmax_coarse,nqmax_fine,nq_blocks_coarse_to_fine, nkmax_coarse, nkmax_fine,nk_blocks_coarse_to_fine,iq_index,nu_index)
        
        command_plus  = build_command(HCK,type_of_integral,mu,T,nkmax_coarse,nkmax_fine, nk_blocks_coarse_to_fine, n_hw,hw_max,Gamma,hw_ph,q_ph1,E_ph1,filename_plus)
        command_minus = build_command(HCK,type_of_integral,mu,T,nkmax_coarse,nkmax_fine, nk_blocks_coarse_to_fine,n_hw,hw_max,Gamma,hw_ph,q_ph2,E_ph2,filename_minus)

        for label, command in zip(['plus','minus'],[command_plus, command_minus]):
            # Take out finished jobs
            while len(set_of_processes) >= max_processes:

                time.sleep(5)
                set_of_finished_processes = set()
                for job in set_of_processes:
                    if job.poll() is not None:
                        set_of_finished_processes.add(job)

                set_of_processes.difference_update(set_of_finished_processes)


            print 'Doing nu =%i, iq_%s = %i of %i q-points pairs'%(nu_index,label,iq_index,number_of_q_pairs )
            print >> log,'Doing nu =%i, iq_%s = %i of %i q-points pairs'%(nu_index,label,iq_index,number_of_q_pairs )
            log.flush()

            job = SUB.Popen(command)
            set_of_processes.add(job)


os.chdir('../')

log.close()
