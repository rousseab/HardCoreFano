import common
reload(common)
from common import *

import scipy.integrate as SI
# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

mpl.rcParams['font.size']='8'




kB = 8.6173324e-5 # eV/K
T  = 300 # K
beta = 1./(kB*T)

mu    =-0.400  # -400 meV
#delta = 0.500  # 500 meV
delta = 0.050  # 50 meV

normalization = density/D_cutoff  

filename ='scattering_spline_mu=%4.3f_eV_delta=%4.3f_eV.nc'%(mu,delta)

SK = ScatteringKernel(mu,beta,delta)
#SK.build_scattering_kernel()
#SK.build_and_write_spline_parameters(filename)
SK.read_spline_parameters(filename)

list_xi = N.arange(-2.*D_cutoff,2*D_cutoff,0.005)

g_plus = get_gamma(list_xi+mu+1j*delta)
g_minus= get_gamma(list_xi+mu-1j*delta)

KR = N.real ( ( g_plus+ g_minus)/2. )
KI = N.real ( (-g_plus+ g_minus)/(2.*1j) )


f  =  function_fermi_occupation(list_xi,0.,beta)
KA =  -2.*N.pi*(hvF)**2/D_cutoff*((list_xi+mu)/D_cutoff)

h = 14  # figure height
w = 22  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(221)
ax2 = fig1.add_subplot(222)
ax3 = fig1.add_subplot(223)
ax4 = fig1.add_subplot(224)

list_ax = [ax1,ax2,ax3,ax4]

ms  = 10 # marker size
mew = 2  # marker edge width
lw  = 3

ax1.plot(list_xi,normalization*KI,'r-',lw=lw,label='$K^I(\\xi+\mu)$')
ax1.plot(list_xi,normalization*KR,'g-',lw=lw,label='$K^R(\\xi+\mu)$')
ax1.plot(list_xi,normalization*KA,'k--',lw=lw,label='$\gamma_A(\\xi+\mu)$')


LI = SK.get_LI_oneD( list_xi )
LR = SK.get_LR_oneD( list_xi )
L   = SK.get_L_oneD( list_xi )

ax2.plot(list_xi,normalization*LI,'r--',lw=lw,label='$L^I(\\xi,0)$')
ax2.plot(list_xi,normalization*LR,'g--',lw=lw,label='$L^R(\\xi,0)$')
ax2.plot(list_xi,normalization*L,'k-',lw=lw,label='$L(\\xi,0)$')

ax2.plot(list_xi,normalization*(f*KA),'k--',lw=lw,label='$L^A(\\xi,0) = f(\\xi) \gamma_A(\\xi+\mu)$')


hw = 0.150 # eV


LI1 = SK.get_LI_oneD( list_xi-hw )
LI2 = SK.get_LI_twoD( list_xi,N.array([hw]) )
dLI = LI1-LI2[:,0]

LR1 = SK.get_LR_oneD( list_xi-hw )
LR2 = SK.get_LR_twoD( list_xi,N.array([hw]) )
dLR = LR1-LR2[:,0]

L1 = SK.get_L_oneD( list_xi-hw )
L2 = SK.get_L_twoD( list_xi,N.array([hw]) )
dL = L1-L2[:,0]

f1  =  function_fermi_occupation(list_xi-hw,0.,beta)
f2  =  function_fermi_occupation(list_xi,0.,beta)
gA  = -2.*N.pi*(hvF)**2/D_cutoff*((list_xi-hw+mu)/D_cutoff)

LA1 = f1*gA
LA2 = f2*gA
dLA =  LA1-LA2

ax3.plot(list_xi,normalization*dLI,'r--',lw=lw,label='$L^I(\\xi-\hbar\omega,0)-L^I(\\xi,\hbar\omega)$')
ax3.plot(list_xi,normalization*dLR,'g--',lw=lw,label='$L^R(\\xi-\hbar\omega,0)-L^R(\\xi,\hbar\omega)$')
ax3.plot(list_xi,normalization*dL,'-',c='gray',alpha=0.5,lw=lw,label='$L(\\xi-\hbar\omega,0)-L(\\xi,\hbar\omega)$')
ax3.plot(list_xi,normalization*dLA,'-',c='black',alpha=1.0,lw=lw,label='$L^A(\\xi-\hbar\omega,0)-L^A(\\xi,\hbar\omega)$')



list_xi2 = N.array([ -mu,0.,mu])
hw = 0.150

L = SK.get_LR_oneD(list_xi)

L2 = SK.get_LR_twoD(list_xi2,N.array([hw,-hw]))

L2_approx = N.array([SK.get_LI_oneD(list_xi2-hw),SK.get_LI_oneD(list_xi2+hw)]).transpose()

for i, xi2 in enumerate(list_xi2): 


    den_plus  = (list_xi-xi2+hw)/((list_xi-xi2+hw)**2+delta**2)
    den_minus = (list_xi-xi2-hw)/((list_xi-xi2-hw)**2+delta**2)

    Y_approx = (L-L2_approx[i,0])*den_plus- (L-L2_approx[i,1])*den_minus

    Y = (L-L2[i,0])*den_plus- (L-L2[i,1])*den_minus

    ax4.plot(list_xi,normalization*Y,'-',alpha=1.0,lw=lw,label='$\\xi_{n_2 \\bf k}$ = %4.2f eV'%xi2)

for ax in list_ax:
    ax.set_xlabel('$\\xi$ (eV)')
    ax.legend(loc=0,fancybox=True,shadow=True)
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)

    ax.legend(loc=0,fancybox=True,shadow=True)


ax3.legend(title='$\hbar\omega$ = %3.1f meV'%(1000*hw), loc=0, ncol=1, fancybox=True,shadow=True,  borderaxespad=0.)

ax4.set_xlabel('$\hbar\omega$ (eV)')



ax1.set_ylabel('$n K(\\xi+\mu)/D$')
ax2.set_ylabel('$n L(\\xi,0)/D$')
ax3.set_ylabel('$n \Delta L(\\xi)/D$')
ax4.set_ylabel('$n Y(\\xi,\hbar\omega)/D$')

# adjust the figure so that there is no wasted space
fig1.suptitle('$\mu$ = %4.3f eV, $\delta$ = %4.3f eV'%(mu,delta))
fig1.subplots_adjust(   left    =   0.07,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()

