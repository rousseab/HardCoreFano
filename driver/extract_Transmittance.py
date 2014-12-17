import sys

module_directory = '/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_5.0/modules/'


sys.path.append(module_directory)

from module_Driver import *

from Scientific.IO.NetCDF import NetCDFFile as Dataset

def RefractionIndexBaF2( list_hw_cm1 ):
        # x is the wavelength, in micro meters
        x = 1e4/list_hw_cm1
        ns = N.sqrt( N.abs(1 + 0.643356*x**2/(x**2-0.057789**2) + 0.506762*x**2/(x**2-0.10968**2) + 3.8261*x**2/(x**2-46.3864**2) ))

        return ns

def build_wedge_to_list(Q_wedge,list_reduced_q):
    """
    The goal of this routine is to identify the index of each wedge k-vector in 
    the list of reduced q vectors.
    """

    tol = 1e-12
    list_wedge_to_list = []

    wedge_to_list = {}

    # Create a large array containing the differences between list_k and list_reduced_q

    x = Q_wedge.list_k[:,0,N.newaxis]-list_reduced_q[N.newaxis,:,0]
    y = Q_wedge.list_k[:,1,N.newaxis]-list_reduced_q[N.newaxis,:,1]

    # the z array contains the distances between k and q vectors
    z = N.sqrt(x**2+y**2)

    # Below, I contains the indices of the k vectors, and J the corresponding indices for the q vectors
    I, J = N.where(z < tol)

    wedge_to_list = {ik:jq for ik,jq in zip(I,J)}


    return wedge_to_list

def phonon_propagator(list_hw,list_hw_ph,width):
    """
    Build the phonon propagator, for

    list_hw_ph ~ [nq]
    list_hw    ~ [nw]

    D(q,w) ~ [nq,nw]
    """

    W = list_hw+1j*width

    Q = list_hw_ph[:,N.newaxis]

    D =  2.*Q/(W**2-Q**2)  

    return D

netcdf_filename = 'HCF_nq=8_32_3_nk=8_32_3_mu=-400_meV_Gamma=50_meV_delta=50_meV.nc'

output_filename = netcdf_filename.replace('nc','dat').replace('HCF','Transmittance')
#================================================================================
# PARAMETERS
#================================================================================

# hyperfine structure constant
alpha = 7.297e-3

# The electron-phonon coupling parameters
u_const = u_a0
v_const = 0.
#v_const = v_a0    

# Phonon width
ph_width = 0.0025 # 2.5 meV

#================================================================================
# Read file and compute
#================================================================================

ncfile   = Dataset(netcdf_filename,'r')

# Extract Q grid parameters
q_nmax_coarse = ncfile.q_nmax_coarse[0]
q_nmax_fine   = ncfile.q_nmax_fine[0]
q_nmax_block  = ncfile.q_nmax_block[0]

mu = ncfile.mu[0]


list_q  = ncfile.variables['list_q'][:]

# Extract other parameters from file
list_hw  = ncfile.variables['list_hw'][:]
list_nu  = ncfile.variables['phonon_nu'][:]
HW_PH    = ncfile.variables['phonon_frequencies'][:,:]

nw = len(list_hw)


# Build the connection between all the q-vectors in wedges and the
# list of reduced q vectors
include_Gamma = True
Q_grid = TesselationDoubleGrid(q_nmax_coarse, q_nmax_fine, q_nmax_block, include_Gamma )

Q_wedge = Q_grid.list_wedges[0]
nk    = len(Q_wedge.list_k)

wedge_to_list = build_wedge_to_list(Q_wedge,list_q)

# We are summing on a single wedge, but there are 12 wedges contributing equally.
wedge_factor = 12.

#================================================================================
# Prepare integration
#================================================================================

PI = complex(0.,0.)*N.zeros(nw)

for inuu, nu in enumerate(list_nu):
    inu = N.int(inuu)
    print 'nu = %i'%nu

    # Extract all data from file for this mode
    # mode frequency for this nu, for every q
    list_hw_ph = HW_PH[inu,:]


    # dimensions: nu, number_of_q_points, number_of_frequencies
    re_HH = ncfile.variables['Re_HH'][inu,:,:]
    im_HH = ncfile.variables['Im_HH'][inu,:,:]

    # introduce the u parameters
    HH = u_const**2*(complex(1.,0.)*re_HH+complex(0.,1.)*im_HH)


    # All arguments of propagator are in eV.
    # the function returns a result in units of 1/eV
    # Convert to fundamental units (Ha)
    Dph   = Ha_to_eV*phonon_propagator(list_hw,list_hw_ph,ph_width)

    DHH   = Dph*HH   
    # The summand 
    Reduced_Summand_PI = -nimp_a02*DHH

    #-------------------------- 
    #
    # PI = 1/N sum_{q} (...)
    #    = Omega/N 1/(2pi)^2 int d2q (...)
    #
    #-------------------------- 
    conversion_factor = Area/(2.*N.pi)**2
    Integral  = 0.

    #TEST_INTEGRAL = 0.

    # perform the integral, to extract the contribution to Pi from this mode


    list_Fk   = complex(0.,0.)*N.zeros([nk,nw])
    for i, k in enumerate(Q_wedge.list_k):
        ir = wedge_to_list[i]
        list_Fk[i,:] = 1.*Reduced_Summand_PI[ir,:] 
    
    Integral += conversion_factor*AreaIntegrator(Q_wedge,list_Fk)


    # Pi in fundamental units
    PI += Integral*wedge_factor  

ncfile.close()
# Compute average conductivity, in atomic units (e^2/hbar)
# Skip first frequency, which is zero

SIGMA = 1j*Ha_to_eV/list_hw[1:]*PI[1:]


ns = RefractionIndexBaF2( list_hw[1:]*eV_to_cm1)
# Compute contribution to transmittance from model 
list_dT = -8.*N.pi*alpha*N.real(SIGMA/(1.+ns))


# Write result to a file

file = open(output_filename ,'w')

print >> file, "#==========================================================="
print >> file, "# Computed value of the relative transmittance of graphene ="
print >> file, "# -------------------------------------------------------- ="
print >> file, "#    Computations performed by Bruno Rousseau              ="
print >> file, "#              CODE: HCF 5.0                               ="
print >> file, "#                                                          ="
print >> file, "#  Formula:                                                ="
print >> file, "#                                                          ="
print >> file, "#  dT(w) = - 8 pi alpha/(1+ns) x                           ="
print >> file, "#                        Re[ sigma_{imp} ]/(e^2/hbar)      ="
print >> file, "#                                                          ="
print >> file, "# Parameters:                                              ="
print >> file, "#                                                          ="
print >> file, "#      mu =  %+4i meV                                      ="%(1000*mu)
print >> file, "#      ns = (BaF2 model, around 1.45)                      ="
print >> file, "#      n_imp x a0^2 = 0.023 * 2/Area  (2.3% impurities)    ="
print >> file, "#==========================================================="
print >> file, "#   hw  (meV)          dT (unitless)                       ="
print >> file, "#==========================================================="

for hw, dT in zip(list_hw[1:], list_dT):
    print >> file, "%12.8f     %+16.12f"%(hw,dT)

file.close()

