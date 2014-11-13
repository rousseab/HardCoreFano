import common
reload(common)
from common import *

#----------------------------------------
# import modules
#----------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties

import matplotlib.cm as cm
import matplotlib.tri as tri


#================================================================================
# Compute set up some parameters
#================================================================================

T = 300. # K
beta = 1./(kB*T)
mu   = -0.400 # eV

Gamma_width = 0.050 # eV

list_hw = N.arange(100.,251.,1.)/1000.

#================================================================================
# prepare the k-grid
#================================================================================
nmax_coarse = 12
nmax_fine   = 48
n_blocks_coarse_to_fine = 3
include_Gamma = False

grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma )

#================================================================================
# Fetch some phonons
#================================================================================
list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()

nu = 5
q_vector = N.array([0.1,0.2])

path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_3.0/modules/'
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)

list_hw_ph, list_hw_renormalized_ph, Eq = Renormalizor.get_frequency_and_polarization(q_vector)

hw_ph = list_hw_renormalized_ph[nu]
E_phonon_polarization = Eq[:,nu]

#================================================================================
# Compute the Loop function
#================================================================================
HQ = Compute_Loop_Function( mu, beta, q_vector, E_phonon_polarization, hw_ph, grid, list_hw, Gamma_width)

HQ.Compute_Hq()

filename = 'test.nc'
write_to_file(HQ,nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,hw_ph,filename)
