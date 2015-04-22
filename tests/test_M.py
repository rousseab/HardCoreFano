#================================================================================
# Plot a nice tesselation of the 1BZ
#
#================================================================================
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

mpl.rcParams['font.size'] = 20.
legendfonts = FontProperties(size=16)


path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_10.0/modules/'
list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)


q_vector = N.random.random(2)

list_hw, list_hw_renormalized, Eq = Renormalizor.get_frequency_and_polarization(q_vector)

nu = N.random.random_integers(0,5)
E_phonon_polarization = Eq[:,nu]
hw_nu_q               = list_hw[nu]

print 'nu = %i'%nu
print 'q  = ( %8.4f,  %8.4f)'%(q_vector[0],q_vector[1])

#nmax = 5
#grid = TesselationGrid(nmax)

nmax_coarse = 8
nmax_fine   = 32
n_blocks_coarse_to_fine = 1
include_Gamma = True

grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma )


for i,wedge in enumerate(grid.list_wedges):

    M_plus  = MGridFunction(q_vector,E_phonon_polarization ,hw_nu_q, wedge)
    M_minus = MGridFunction(-q_vector,N.conjugate(E_phonon_polarization), hw_nu_q, wedge)


    if i == 0:
        list_k = wedge.list_k

        list_M_plus = N.transpose(M_plus.M,(2,0,1))
        list_M_minus= N.transpose(M_minus.M,(2,0,1))

        list_J_x_plus = N.transpose(M_plus.J_x,(2,0,1))
        list_J_x_minus= N.transpose(M_minus.J_x,(2,0,1))

        list_J_y_plus = N.transpose(M_plus.J_y,(2,0,1))
        list_J_y_minus= N.transpose(M_minus.J_y,(2,0,1))


    else:
        list_k = N.concatenate([list_k,wedge.list_k])

        list_M_plus = N.concatenate([list_M_plus ,N.transpose(M_plus.M,(2,0,1))])
        list_M_minus= N.concatenate([list_M_minus,N.transpose(M_minus.M,(2,0,1))])

        list_J_x_plus = N.concatenate([list_J_x_plus ,N.transpose(M_plus.J_x,(2,0,1))])
        list_J_x_minus= N.concatenate([list_J_x_minus,N.transpose(M_minus.J_x,(2,0,1))])

        list_J_y_plus = N.concatenate([list_J_y_plus ,N.transpose(M_plus.J_y,(2,0,1))])
        list_J_y_minus= N.concatenate([list_J_y_minus,N.transpose(M_minus.J_y,(2,0,1))])


# Find the k and minus k

list_i_plus  = []
list_i_minus = []

tol = 1e-12
for i, kplus in  enumerate(list_k):
    list_i_plus.append(i)

    for j, kminus in enumerate(list_k):
        if N.linalg.norm(kplus+kminus) < tol:
            list_i_minus.append(j)
            # you have to break; the k points are 
            # repeated! more than one solutions will be found
            break

found_Error = False
for ip, im in zip(list_i_plus,list_i_minus):

    Mp = list_M_plus[ip]
    Mm = list_M_minus[im]

    Jxp = list_J_x_plus[ip]
    Jxm = list_J_x_minus[im]

    Jyp = list_J_y_plus[ip]
    Jym = list_J_y_minus[im]


    error_M = N.linalg.norm( Mp+N.conjugate(Mm) )
    error_Jx= N.linalg.norm( Jxp+N.conjugate(Jxm) )
    error_Jy= N.linalg.norm( Jyp+N.conjugate(Jym) )

    if error_M > 1e-10:
        print 'ERROR IN M SYMMETRY! error = %12.6e'%error
        found_Error = True

    if error_Jx> 1e-10:
        print 'ERROR IN Jx SYMMETRY! error = %12.6e'%error_Jx
        found_Error = True

    if error_Jy> 1e-10:
        print 'ERROR IN Jy SYMMETRY! error = %12.6e'%error_Jy
        found_Error = True



if not found_Error:
    print 'symmetries are respected'


