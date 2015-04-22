#================================================================================
# Plot a nice tesselation of the 1BZ
#
#================================================================================
import common
reload(common)
from common import *


import matplotlib.pyplot as plt
import seaborn; seaborn.set(font_scale=2)

path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_10.0/modules/'
list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)

q_vector = N.array([0.,0.])

list_hw, list_hw_renormalized, Eq = Renormalizor.get_frequency_and_polarization(q_vector)


print 'q  = ( %8.4f,  %8.4f)'%(q_vector[0],q_vector[1])


print "# ====  Computing  ===="
print "#"
print "#  K^{alpha,nu}(omega) = 1/Omega sum_{k sigma} J^{alpha}_{n1n2}(k) g^{nu}_{n2n1}(k,k)  f(xi_{n1 k}) - f(xi_{n2 k})"
print "#                                                                                      ---------------------------"
print "#                                                                                      hw +i delta+xi_{n1 k}-xi_{n2 k}"
print "#====================="


"""
nmax_coarse = 32
nmax_fine   = 128
n_blocks_coarse_to_fine = 12
"""
nmax_coarse = 16
nmax_fine   = 64
n_blocks_coarse_to_fine = 6


clip_grid = False
include_Gamma = False
grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma , \
                                clip_grid=clip_grid,clip_energy=0.)

delta   = 0.250 # eV
mu      = -0.500 # eV
T       = 300  # K
beta    = 1./(kB*T)

list_hw_ext =  N.arange(-2.,2.,0.01)

integration_factor = 1./(2.*N.pi**2)

K = complex(0.,0.)*N.zeros([2,6,len(list_hw_ext)])

for nu in N.arange(0,6):
    print 'nu = %i'%nu
    E_phonon_polarization = Eq[:,nu]
    hw_nu_q               = list_hw[nu]

    for i, wedge in enumerate(grid.list_wedges):

        M = MGridFunction(q_vector,E_phonon_polarization ,hw_nu_q, wedge)

        for n1 in [0,1]:
            eps_n1k = M.epsilon_k[n1,:]
            f_n1k   = function_fermi_occupation(eps_n1k,mu,beta)
            for n2 in [0,1]:
                eps_n2k = M.epsilon_k[n2,:]
                f_n2k   = function_fermi_occupation(eps_n2k,mu,beta)

                CHI =  (f_n1k[:,N.newaxis]-f_n2k[:,N.newaxis])/(eps_n1k[:,N.newaxis]-eps_n2k[:,N.newaxis]+1j*delta+list_hw_ext[N.newaxis,:])

                #print 'n1 = %i, n2 = %i MAX[CHI] = %10.4f'%(n1,n2,N.abs(CHI).max())

                gnu = u_a0*N.abs(M.g_nu_u[n2,n1,:]) # in funamental units, Ha
                #gnu = M.g_nu_u[n2,n1,:]

                for i_alpha, J_alpha in zip([0,1],[M.J_x,M.J_y]):
                    J = J_alpha[n1,n2,:] # in fundamental units,  e hbar/ m a0

                    gJ = (gnu*Ha_to_eV)*J

                    list_Fk = gJ[:,N.newaxis]*CHI
                    #K[i_alpha,nu,:] += N.dot(CHI,gnu*J)

                    # in units of e hbar /m a0^3
                    K[i_alpha,nu,:] += integration_factor*AreaIntegrator(wedge,list_Fk)



w = 14
h = 8
fig1 = plt.figure(figsize=(w,h))
ax1  = fig1.add_subplot(121)
ax2  = fig1.add_subplot(122)

for ax in [ax1,ax2]:
    ax.set_xlabel(' $\hbar\omega$ (eV)')

ax1.set_ylabel(' Re[K] ($e \hbar /m a_0^3$)')
ax2.set_ylabel(' Im[K] ($e \hbar /m a_0^3$)')

nu = 5

for i_alpha in [0,1]:
    Y = K[i_alpha,nu,:]
    ax1.plot(list_hw_ext,N.real(Y),label='i_alpha = %i, nu = %i'%(i_alpha,nu))
    ax2.plot(list_hw_ext,N.imag(Y),label='i_alpha = %i, nu = %i'%(i_alpha,nu))


ax1.legend(loc=0)
plt.tight_layout()
plt.show()

