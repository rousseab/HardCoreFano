#================================================================================
#
#           Module Constants
#           ===========================
#
#   This "module" will simply contain, once and for all, various geometric
#   constants relating to Graphene. It also imports all the necessary
#   modules.
#================================================================================

import numpy as N
import scipy as S
import scipy.special as SS
from copy import deepcopy


# Make sure we have no monkey business
#import warnings
#warnings.simplefilter("error", RuntimeWarning)

#----------------------------------------
# universal constants 
#----------------------------------------

bohr_in_angst  = 0.529177249   # conversion between angstrom and bohr
Ha_to_eV       = 27.211396132  # conversion between Hartree and eV
eV_to_cm1      = 8065.73
kB             = 8.6173324e-5  # Boltzmann constant, in eV / K

#----------------------------------------
# parameters specific to graphene
#----------------------------------------
acell_angstrom = 2.461                  # in Angstrom
acell = acell_angstrom/bohr_in_angst    # in bohrs

twopia = 2.*N.pi/acell
Area   = N.sqrt(3.)/2.*acell**2

# latice vectors
a1 = acell*N.array([1.,0.]) 
a2 = acell*N.array([0.5,N.sqrt(3)/2.]) 

# reciprocal lattice vectors
b1 = 2.*N.pi/acell*N.array([1.,-N.sqrt(3.)/3.]) 
b2 = 2.*N.pi/acell*N.array([0.,2.*N.sqrt(3.)/3.])

# basis vectors, describing carbon positions, in units of a
tau_A = N.array([1./2., N.sqrt(3.)/6.])*acell
tau_B = N.array([1./2.,-N.sqrt(3.)/6.])*acell

# vectors connecting near neighbors in graphene
tau1  = 0.5*acell*N.array([1.,N.sqrt(3.)/3.])
tau2  = 0.5*acell*N.array([-1.,N.sqrt(3.)/3.])
tau3  = 0.5*acell*N.array([0.,-2.*N.sqrt(3.)/3.])

tau1_minus_tau3 = tau1 - tau3  
tau2_minus_tau3 = tau2 - tau3  

hat_tau1 = N.array([ N.sqrt(3.)/2.,0.5])
hat_tau2 = N.array([-N.sqrt(3.)/2.,0.5])
hat_tau3 = N.array([ 0.0,-1.0])

norm_K_point = 4./3.*N.pi/acell
norm_M_point = N.sqrt(3.)/2.*norm_K_point 

# tight binding parameter 
tight_binding_gamma1 =  -2.7 # eV; negative according to Peres et al.
hvF = N.sqrt(3.)/2.*acell*N.abs(tight_binding_gamma1)  # hbar times Fermi velocity, in bohr x eV

# Estimating the various unit carrying quantities

proton_mass   = 1.67262178e-27 # kg
electron_mass = 9.10938291e-31 # kg
carbon_mass   = 12.0107*proton_mass

impurity_density = 0.023 # percentage of sites with defects, N_imp / N_carbon, 2.3%

density  = 2./Area # two carbon atoms per unit cell area
nimp_a02 = impurity_density*density  

me_M     = electron_mass/carbon_mass

# Values fitted from the literature
u_a0_Gamma =  2.80
u_a0_K     =  3.47
u_a0    = ( u_a0_Gamma + u_a0_K)/2. # eV; average value from GW results of Mauri et al.
v_a0    = 0.0   # eV

