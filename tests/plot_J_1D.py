#================================================================================
#
# This test computes the integral for J numerically, and plots the 
# results for various values of the parameters xi_1, xi_3, hw.
#
#================================================================================
import common
reload(common)
from common import *

import module_numerical_tests 
reload(module_numerical_tests )
from module_numerical_tests import *

import os
import numpy as N
import matplotlib.pyplot as plt
import seaborn; seaborn.set(font_scale=2)


#================================================================================
#parameters 
#================================================================================
T = 300.
beta = 1./(kB*T)

kernel_Gamma_width = 0.200 # eV
Green_Gamma_width  = 0.100
mu = -0.400
hw_ext = 0.150


#state_n2 = 'interband'
state_n2 = 'intraband'

# read numerical integration result
filename = 'J_Green_Gamma=%i_meV_mu=%i_meV_hw=%3.1f_meV.nc'%(1000*Green_Gamma_width,1000*mu,1000*hw_ext)
handler = NetcdfHandler(filename)
handler.open_ncfile(mode='r')
list_xi_1, list_xi_3, mu, Green_Gamma_width, kernel_Gamma_width, hw_ph, beta = handler.read_attributes()
J_numerical = handler.read_J_array('numerical_J_%s'%state_n2)
handler.close_ncfile()


i_xi_1 = 100

xi_1 = list_xi_1[i_xi_1]
list_J_numerical = J_numerical[i_xi_1,:] 

parameters = {      '$\Gamma_\gamma$': kernel_Gamma_width,
                    '$\Gamma_G$': Green_Gamma_width ,
                    '$\mu$': mu,
                    '$\hbar\omega$': hw_ext,
                    '$\\xi_{n_1 k}$': xi_1}



# Compute along a given line

spline_filename = 'scattering_spline_mu=%i_meV_Green_Gamma_width=%i_meV.nc'%(1000*mu,1000*Green_Gamma_width)
os.system('cp %s scattering_spline.nc'%spline_filename )

q_vector = N.array([0.,0.])
wedge    = DummyWedge(N.array([q_vector])) # just a hack dummy argument
Jobject  = JGridFunction(q_vector, hw_ext, Green_Gamma_width, kernel_Gamma_width, mu, beta, wedge)


list_J_code = complex(0.,0.)*N.zeros_like(list_J_numerical)

tol = 1e-8
delta = 1e-4
xi_1_array = N.array([xi_1])

if state_n2 == 'intraband':
    xi_2_array = xi_1_array 
elif state_n2 == 'interband':
    xi_2_array = -xi_1_array -2*mu

for i3, xi_3 in enumerate(list_xi_3):
    print 'doing i3 = %i'%i3
    xi_3_array = N.array([xi_3])

    if N.abs(xi_1 - xi_3) < tol:
        D_1 = Jobject.get_D(xi_1_array, xi_2_array)[0]
        D_3 = Jobject.get_D(xi_3_array+delta, xi_2_array)[0]

        list_J_code[i3] += (D_1-D_3)/(-delta)
    else:

        D_1 = Jobject.get_D(xi_1_array, xi_2_array)[0]
        D_3 = Jobject.get_D(xi_3_array, xi_2_array)[0]
        
        list_J_code[i3] += (D_1-D_3)/(xi_1-xi_3)





h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(221)
ax2 = fig1.add_subplot(222)

ax3 = fig1.add_subplot(223)
ax4 = fig1.add_subplot(224)


kwargs = { 'ms':8,'lw':4}

dJ = list_J_numerical-list_J_code

ax1.plot(list_xi_3, N.real(list_J_numerical),'-',label='Numerical',**kwargs)
ax1.plot(list_xi_3, N.real(list_J_code),'o',label='Code',**kwargs)

ax2.plot(list_xi_3, N.imag(list_J_numerical),'-',label='Numerical',**kwargs)
ax2.plot(list_xi_3, N.imag(list_J_code),'o',label='Code',**kwargs)

ax3.plot(list_xi_3, N.real(dJ),'-',**kwargs)
ax4.plot(list_xi_3, N.imag(dJ),'-',**kwargs)



x1 = [xi_1,xi_1]
x2 = [xi_1+hw_ext,xi_1+hw_ext]
x3 = [xi_1-hw_ext,xi_1-hw_ext]


ax1.legend(loc=0)
for ax in [ax1,ax2,ax3,ax4]:
    ax.set_xlabel('$\\xi_3$ (eV)')
    ax.set_xlim([-1.,1.])


    ymin, ymax = ax.set_ylim()

    for x in [x1,x2,x3]:
        ax.plot(x,[ymin,ymax],'k--',label='__nolabel__')
    


ax1.set_ylabel('Re[$J$]')
ax2.set_ylabel('Im[$J$]')

ax3.set_ylabel('Re[$J_{num}-J_{code}$]')
ax4.set_ylabel('Im[$J_{num}-J_{code}$]')


title = ''
for str,e in parameters.iteritems():
    title += '%s = %3.1f meV  '%(str,1000*e)

plt.suptitle(title)


fig1.subplots_adjust(   left    =       0.10,
                        bottom  =       0.10,
                        right   =       0.90,
                        top     =       0.90,
                        wspace  =       0.30,
                        hspace  =       0.30)

#plt.tight_layout()

plt.show()
