import common
reload(common)
from common import *


# parameters

kB = 8.6173324e-5 # eV/K
T  = 300 # K
beta = 1./(kB*T)

mu    = -0.400 # meV 
kernel_Gamma_width = 0.200  # 200 meV
Green_Gamma_width  = 0.025   

list_hw_ext = N.arange(0.000,0.220,0.001)

filename ='scattering_spline_mu=%i_meV_Green_Gamma_width=%i_meV.nc'%(1000*mu,1000*Green_Gamma_width) 
n_xi_grid = 8001

SK = ScatteringKernel(mu, beta, kernel_Gamma_width, Green_Gamma_width, spline_order = 1)
SK.build_scattering_kernel_integrals(list_hw_ext, n_xi_grid = n_xi_grid )
SK.build_and_write_spline_parameters(filename)
#SK.read_spline_parameters(filename)

