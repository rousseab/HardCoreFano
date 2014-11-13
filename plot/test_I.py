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


q_vector = N.random.random(2)

list_hw  = N.arange(0.,300.,1.)/1000. # in eV

T = 300. # K

beta = 1./(kB*T)
mu   = -0.400 # eV


delta_width = 0.050 # eV
nmax = 6
grid = TesselationGrid(nmax)

for i,wedge in enumerate(grid.list_wedges):
       	
	I_plus  = IGridFunction( q_vector,list_hw, delta_width, mu, beta, wedge)
	I_minus = IGridFunction(-q_vector,list_hw, delta_width, mu, beta, wedge)

	if i == 0:
		list_k = wedge.list_k

		list_I_plus = N.transpose(I_plus.I,(1,0,2))
		list_I_minus= N.transpose(I_minus.I,(1,0,2))

	else:
		list_k = N.concatenate([list_k,wedge.list_k])

		list_I_plus = N.concatenate([list_I_plus ,N.transpose(I_plus.I,(1,0,2))])
		list_I_minus= N.concatenate([list_I_minus,N.transpose(I_minus.I,(1,0,2))])

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

	Ip = list_I_plus[ip]
	Im = list_I_minus[im]

	error = N.linalg.norm( Ip-Im) 
	if error > 1e-10: 
		print 'ERROR IN I SYMMETRY! error = %8.4e'%error

		found_Error = True


if not found_Error:
	print 'symmetries are respected'


