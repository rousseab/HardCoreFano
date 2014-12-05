# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.misc import pade


# Define some fake data, for the sake of this example
acell_angstrom = 2.461                  # in Angstrom
bohr_in_angst = 0.5291772083
acell = acell_angstrom/bohr_in_angst    # in bohrs
tight_binding_gamma1 =  -2.7 # eV; negative according to Peres et al.
hvF = N.sqrt(3.)/2.*acell*N.abs(tight_binding_gamma1)  # hbar times Fermi velocity, in bohr x eV 

twopia = 2.*N.pi/acell
Area   = N.sqrt(3.)/2.*acell**2
kc = N.sqrt(2.*N.pi/Area)
D = 2*hvF/acell*(N.pi**2/3.)**0.25

print 'D = %3.1f eV'%D

def gamma(z):
    
    prefactor = 2*N.pi*(hvF)**2/D
    #prefactor = 1.

    x = D/z
    g = prefactor*x/N.log(1.-x**2)
    return g

def gamma_inf(z,nmax=1):
    
    prefactor = -2*N.pi*(hvF)**2/D
    #prefactor = -1

    u = z/D

    series = complex(0.,0.)*N.zeros_like(u)
            
    for n in N.arange(nmax+1):
        pow = 2*(nmax-n)
        series += u**pow/(n+1)

    g = prefactor*u**(2*nmax+1)/series

    return g

def delta_gamma_model(z):
    
    prefactor = -2*N.pi*(hvF)**2/D
    #prefactor = -1.

    u = z/D

    g = -prefactor/u/2.

    return g

kB = 8.6173324e-5 # eV/K
T  = 300 # K

beta = 1./(kB*T)

mu = -0.4
list_delta = [-0.05,0.025,0.05]
x = N.arange(-0.8*D,0.8*D,0.01)


h = 10  # figure height
w = 10  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(122)
ax2 = fig1.add_subplot(121)

ms  = 10 # marker size
mew = 2  # marker edge width


ms  = 8
mew = 2
lw  = 4

for delta,ls in zip(list_delta,['-','-','--']):
    g = gamma(x+1j*delta)
    ni_sigma = 2./Area*g/D
    im = N.imag(ni_sigma)
    re = N.real(ni_sigma)


    ax1.plot(x/D,re/D,ls=ls,lw=lw)
    ax2.plot(x/D,-im/D,ls=ls,lw=lw,label='$\delta$ = %4.1f meV'%(1000*delta))



ax1.plot([mu/D,mu/D],[-3,3],'k--',lw=lw)
ax2.plot([mu/D,mu/D],[0,10],'k--',lw=lw,label='$\mu$ = %3.1f meV'%(1000*mu))



for ax in [ax1,ax2]:
    ax.set_xlabel('Re[x]/D')
    ax.legend(loc=0,fancybox=True,shadow=True)
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)

ax1.set_ylabel('Real[$\Sigma_{imp}$]/(D $N_{imp}/N$)')
ax2.set_ylabel('-Imag[$\Sigma_{imp}$]/(D $N_{imp}/N$)')


ax1.set_ylim([-3.,3.])
ax2.set_ylim([-10.,10.])
fig1.suptitle('Scattering kernel function: Real axis')

# adjust the figure so that there is no wasted space
plt.subplots_adjust(    left    =   0.10,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()
