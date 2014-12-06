import common
reload(common)
from common import *
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm



T    = 300.
beta = 1./(kB*T)

mu = -0.400 # eV

Gamma = 0.0050 # eV

K1 = 2./3.*twopia*N.array([1.,0.])

path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_4.0/modules/'

list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)
q_vector = K1+N.random.random(2)/50*twopia

list_hw, list_hw_renormalized, Eq = Renormalizor.get_frequency_and_polarization(q_vector)

nu = 5
E_phonon_polarization = Eq[:,nu]
hw_nu_q               = list_hw[nu]

print 'nu = %i'%nu
print 'q  = ( %8.4f,  %8.4f)'%(q_vector[0],q_vector[1])

list_hw = N.arange(0,0.250,0.001) # eV


IHq_plus = Compute_Imaginary_Loop_Function(mu, beta, q_vector, E_phonon_polarization, hw_nu_q, list_hw,Gamma)
IHq_plus.Compute_Hq() 
IHq_minus= Compute_Imaginary_Loop_Function(mu, beta,-q_vector, N.conjugate(E_phonon_polarization), hw_nu_q, list_hw,Gamma)
IHq_minus.Compute_Hq() 

Hq_minus = IHq_minus.Hq
Hq_plus  = IHq_plus.Hq

# Test symmetry
error = N.linalg.norm( N.conjugate(Hq_plus) - Hq_minus  )

print '|H(q,w)| = %12.8e   error: %12.8e'%(N.linalg.norm(Hq_plus),error)

#list_K_plus, list_K_minus = IHq.build_list_K(I_solutions,list_cos_theta,list_sin_theta)
   
h =  8  # figure height
w = 12  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)
   
HW = N.arange(-0.25,0.25,0.001)

difference = N.zeros_like(HW)

KP = []
KM = []
for ax, IHq, sign in zip([ax1,ax2],[IHq_plus,IHq_minus],[1.,-1]):

    I_solutions, I_solutions_upper_band, \
    I_solutions_lower_band, list_cos_theta, list_sin_theta = IHq.get_cos_theta(HW)
    F  = IHq.get_frequency_prefactor(HW,I_solutions,list_sin_theta)

    list_K_plus, list_K_minus = IHq.build_list_K(I_solutions,list_cos_theta,list_sin_theta)

    KP.append(list_K_plus)
    KM.append(list_K_minus)

    difference += sign*F

    eQ = IHq_plus.nQ*hvF
    ax.plot(HW,F,'k--',lw=2,label='all')
    ax.plot(HW[I_solutions_lower_band],F[I_solutions_lower_band]**2,'r-',lw=4,label='lower band')
    ax.plot(HW[I_solutions_upper_band],F[I_solutions_upper_band]**2,'g-',lw=4,label='upper band')

    y = 60000
    x = 2*mu+eQ
    ax.plot([x,x],[-y,y],'--',c='grey',lw=4,label='expected divergences')
    x = 2*mu-eQ
    ax.plot([x,x],[-y,y],'--',c='grey',lw=4,label='__nolabel__')

    x = -eQ
    ax.plot([x,x],[-y,y],'--',c='grey',lw=4,label='__nolabel__')

    x =  eQ
    ax.plot([x,x],[-y,y],'--',c='grey',lw=4,label='__nolabel__')


    ax.set_xlim([-0.25,0.25])
    ax.legend(loc=0)
plt.show()

