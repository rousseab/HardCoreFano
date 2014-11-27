import common
reload(common)
from common import *

import scipy.integrate as SI
# import the usual modules
import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm



kB = 8.6173324e-5 # eV/K
T  = 300 # K
beta = 1./(kB*T)

mu =   -0.40 # -400 meV

normalization = density/D_cutoff  


delta = 0.010 # 50 meV

n_u   = 4001
u_max = 4.*D_cutoff
u_min =-4.*D_cutoff
d_u   = (u_max-u_min)/(n_u-1.)
iloop  = N.arange(n_u)
list_u = u_min+d_u*iloop

KK = KramersKronig(list_u)



h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax3 = fig1.add_subplot(133)

list_ax = [ax1,ax2,ax3]

ms  = 10 # marker size
mew = 2  # marker edge width
lw  = 3


gC_plus    = get_gammaC(list_u+mu+1j*delta)
gC_minus   = get_gammaC(list_u+mu-1j*delta)
KI = (gC_minus-gC_plus)/(2.*1j)
KR = (gC_minus+gC_plus)/2.


list_KK_KI = N.real(1j*KK.apply_kramers_kronig_FFT_convolution(KI))


ax1.plot(list_u,normalization*N.real(KI),'b-',lw=lw,label='$K_C^I(\epsilon+\mu)$')
ax1.plot(list_u,normalization*N.real(KR),'r-',lw=lw,label='$K_C^R(\epsilon+\mu)$')
ax1.plot(list_u,-normalization*N.real(list_KK_KI),'g--',lw=lw,label='$-KK(K_C^I(\epsilon+\mu))$')


f  = function_fermi_occupation(list_u,0.,beta)

list_FKR = f*KR
list_Fi = f*KI 
list_Fr = N.real(1j*KK.apply_kramers_kronig_FFT_convolution(list_Fi))

ax2.plot(list_u,normalization*list_Fi,'b-',lw=lw,label='$f(\epsilon) K_C^I(\mu+\epsilon)$')
ax2.plot(list_u,normalization*list_Fr,'r-',lw=lw,label='$L(f(\epsilon) K_C^I(\mu+\epsilon))$')


ax3.plot(list_u,normalization*list_Fr,'r-',lw=lw,label='$L(f(\epsilon) K_C^I(\mu+\epsilon))$')
ax3.plot(list_u,normalization*list_FKR,'b-',lw=lw,label='$f(\epsilon) K_C^R(\mu+\epsilon)$')
ax3.plot(list_u,normalization*(list_Fr+list_FKR),'g--',lw=lw,label='sum')


SK = ScatteringKernel(mu,beta,delta)

filename ='scattering_spline.nc'
#SK.build_scattering_kernel()
#SK.build_and_write_spline_parameters(filename)
SK.read_spline_parameters(filename)

list_xi = N.arange(-4*D_cutoff,4*D_cutoff,0.01)
list_SK = SK.get_LC(list_xi)
list_SK2 = SK.get_L(list_xi)

#ax3.plot(SK.list_x,normalization*(SK.LC_kernel),'k--',lw=lw,label='CHECK')
ax3.plot(list_xi,normalization*(list_SK),'k--',lw=lw,label='LC')
ax3.plot(list_xi,normalization*(list_SK2),'b--',lw=lw,label='L')


def integrand(y,x):

    fy = function_fermi_occupation(N.array([y]),0.,beta)

    gC_plus    = get_gammaC(y+mu+1j*delta)
    gC_minus   = get_gammaC(y+mu-1j*delta)

    KI = (gC_minus-gC_plus)/(2.*1j)

    P  = (y-x)/((y-x)**2+delta**2)

    integrand = P/N.pi*KI*fy

    return N.real(integrand[0])

"""
list_u_small      = []
list_Fr_numerical = []
list_err = []
for u in x[::40]:

    list_u_small.append(u)
    im = SI.quad(integrand, -1.5*D_cutoff, 1.5*D_cutoff, args=(u))
    list_Fr_numerical.append(im[0])
    list_err.append(im[1])


list_u_small = N.array(list_u_small)
list_Fr_numerical = N.array(list_Fr_numerical)
ax2.plot(list_u_small,normalization*list_Fr_numerical ,'g--',lw=lw,label='$L(f K_C^I)$')
"""



for ax in list_ax:
    ax.set_xlabel('$\epsilon$ (eV)')
    ax.legend(loc=0,fancybox=True,shadow=True)
    ax.grid(True,linestyle='-',color='grey',alpha=0.5)

    ax.legend(loc=0,fancybox=True,shadow=True)


ax1.set_ylabel('$n K(\epsilon)/D$')

# adjust the figure so that there is no wasted space
fig1.suptitle('$\mu$ = %4.3f eV, $\delta$ = %4.3f eV'%(mu,delta))
fig1.subplots_adjust(   left    =   0.07,
                        bottom  =   0.15,
                        right   =   0.95,
                        top     =   0.95,
                        wspace  =   None, 
                        hspace  =   None)

plt.show()
