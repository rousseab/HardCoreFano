#================================================================================
#
# This test computes the integral for J numerically, and plots the 
# results for various values of the parameters xi_1, xi_3, hw.
#
#================================================================================

import common
reload(common)
from common import *
import os

from module_numerical_tests import *
# import the usual modules
import numpy as N
import matplotlib.pyplot as plt

"""
Green_Gamma_width = 0.100 # eV
mu = -0.400 # eV
hw_ext = 0.150 # eV

Green_Gamma_width = 0.025 # meV
mu = -0.400 # meV
hw_ext = 0.1755 # meV
"""
Green_Gamma_width = 0.025 # meV
mu = -0.400 # meV
hw_ext = 0.200 # meV





#state_n2 = 'interband'
state_n2 = 'intraband'

filename = 'J_Green_Gamma=%i_meV_mu=%i_meV_hw=%3.1f_meV.nc'%(1000*Green_Gamma_width,1000*mu,1000*hw_ext)
handler = NetcdfHandler(filename)
handler.open_ncfile(mode='r')
list_xi_1, list_xi_3, mu, Green_Gamma_width, kernel_Gamma_width, hw_ph, beta = handler.read_attributes()
J_numerical = handler.read_J_array('numerical_J_%s'%state_n2)
handler.close_ncfile()


spline_filename = 'scattering_spline_mu=%i_meV_Green_Gamma_width=%i_meV.nc'%(1000*mu,1000*Green_Gamma_width)
os.system('cp %s scattering_spline.nc'%spline_filename )

q_vector = N.array([0.,0.])
wedge    = DummyWedge(N.array([q_vector])) # just a hack dummy argument
Jobject  = JGridFunction(q_vector, hw_ph, Green_Gamma_width, kernel_Gamma_width, mu, beta, wedge)

J_code = complex(0.,0.)*N.zeros_like(J_numerical)


tol = 1e-8
delta = 1e-4
for i1, xi_1 in enumerate(list_xi_1):
    print 'i1 = %i'%i1
    xi_1_array = N.array([xi_1])

    if state_n2 == 'intraband':
        xi_2_array = xi_1_array 
    elif state_n2 == 'interband':
        xi_2_array = -xi_1_array -2*mu
    
    for i3, xi_3 in enumerate(list_xi_3):
        xi_3_array = N.array([xi_3])

        if N.abs(xi_1 - xi_3) < tol:
            D_1 = Jobject.get_D(xi_1_array, xi_2_array)[0]
            D_3 = Jobject.get_D(xi_3_array+delta, xi_2_array)[0]

            J_code[i1,i3] += (D_1-D_3)/(-delta)
        else:

            D_1 = Jobject.get_D(xi_1_array, xi_2_array)[0]
            D_3 = Jobject.get_D(xi_3_array, xi_2_array)[0]
            
            J_code[i1,i3] += (D_1-D_3)/(xi_1-xi_3)


h = 10  # figure height
w = 18  # figure width

fig = plt.figure(figsize=(w,h))
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(234)
ax3 = fig.add_subplot(232)
ax4 = fig.add_subplot(235)
ax5 = fig.add_subplot(233)
ax6 = fig.add_subplot(236)


L = N.concatenate([J_numerical.flatten(),J_code.flatten()])

vmin = N.min([N.real(L).min(),N.imag(L).min()])
vmax = N.max([N.real(L).max(),N.imag(L).max()])

extent= [list_xi_3.min(),list_xi_3.max(),list_xi_1.min(),list_xi_1.max()]


ax1.set_title('exact J: Real part')
ax2.set_title('exact J: Imaginary part')

ax3.set_title('J CODE: Real part')
ax4.set_title('J CODE: Imaginary part')


ax5.set_title('DIFFERENCE: Real part')
ax6.set_title('DIFFERENCE: Imaginary part')

dJ = J_numerical-J_code

list_ax = [ax1,ax2,ax3,ax4,ax5,ax6]
for ax in list_ax:
    if state_n2 == 'intraband':
        ax.plot(list_xi_1, list_xi_1,'k--')
        ax.plot(list_xi_1, list_xi_1+hw_ext,'k--')
        ax.plot(list_xi_1, list_xi_1-hw_ext,'k--')
    elif state_n2 == 'interband':
        ax.plot(list_xi_1, -list_xi_1-2*mu,'k--')
        ax.plot(list_xi_1, -list_xi_1-2*mu+hw_ext,'k--')
        ax.plot(list_xi_1, -list_xi_1-2*mu-hw_ext,'k--')


    ax.set_xlabel('$\\xi_3$ (eV)')
    ax.set_ylabel('$\\xi_1$ (eV)')


im1 = ax1.imshow(N.real(J_numerical[::-1,:]),vmin=vmin,vmax=vmax,aspect='equal', extent=extent)
im2 = ax2.imshow(N.imag(J_numerical[::-1,:]),vmin=vmin,vmax=vmax,aspect='equal', extent=extent)
im3 = ax3.imshow(N.real(J_code[::-1,:]),vmin=vmin,vmax=vmax,aspect='equal', extent=extent)
im4 = ax4.imshow(N.imag(J_code[::-1,:]),vmin=vmin,vmax=vmax,aspect='equal', extent=extent)
im5 = ax5.imshow(N.real(dJ[::-1,:]),aspect='equal', extent=extent)
im6 = ax6.imshow(N.imag(dJ[::-1,:]),aspect='equal', extent=extent)

list_im = [im1,im2,im3,im4,im5,im6]

for ax,im in zip(list_ax,list_im):
    fig.colorbar(im, ax=ax,shrink=0.75)

plt.tight_layout()
#plt.show()
plt.savefig('test_Green_Gamma=%i_meV_mu=%i_meV_hw=%3.1f_meV_%s.pdf'%(1000*Green_Gamma_width,1000*mu,1000*hw_ext,state_n2))
