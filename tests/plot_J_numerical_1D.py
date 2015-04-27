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

"""
mu = -0.400
Green_Gamma_width  = 0.100
hw_ext = 0.150
"""
Green_Gamma_width  = 0.0250
mu = -0.400
hw_ext = 0.175


Integrand_generator = NumericalJ( hw_ext, mu, beta, kernel_Gamma_width, Green_Gamma_width)

global_integrand_J_f  = Integrand_generator.get_exact_J_f_integrand
global_integrand_J_df = Integrand_generator.get_exact_J_df_integrand

state_n2 = 'intraband'
#state_n2 = 'interband'

xi_1 = 0.0

list_xi_3 = N.arange(-1,1.01,0.05)


parameters = {      '$\Gamma_\gamma$': kernel_Gamma_width,
                    '$\Gamma_G$': Green_Gamma_width ,
                    '$\mu$': mu,
                    '$\hbar\omega$': hw_ext,
                    '$\\xi_{n_1 k}$': xi_1}


tol = 1e-8
delta = 1e-4

if state_n2 == 'intraband':
    xi_2 = xi_1
elif state_n2 == 'interband':
    xi_2 = -xi_1 - 2*mu


global_parameters = {'xi_n1k':xi_1,'xi_n2k':xi_2}

list_re_J_f, list_im_J_f   = Drive_Numerical_Integrals(list_xi_3, global_integrand_J_f, 'xi_n3kq', global_parameters)
list_re_J_df, list_im_J_df = Drive_Numerical_Integrals(list_xi_3, global_integrand_J_df, 'xi_n3kq', global_parameters)



h = 10  # figure height
w = 16  # figure width

fig1 = plt.figure(figsize=(w,h))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)


kwargs = { 'ms':8,'lw':4}

list_re_J = list_re_J_df+list_re_J_f
list_im_J = list_im_J_df+list_im_J_f

ax1.plot(list_xi_3, list_re_J_f,'o-',label='f',**kwargs)
ax1.plot(list_xi_3, list_re_J_df,'o-',label='df',**kwargs)
ax1.plot(list_xi_3, list_re_J,'h-',label='sum',**kwargs)


ax2.plot(list_xi_3, list_im_J_f,'o-',label='f',**kwargs)
ax2.plot(list_xi_3, list_im_J_df,'o-',label='df',**kwargs)
ax2.plot(list_xi_3, list_im_J,'h-',label='sum',**kwargs)

x1 = [xi_1,xi_1]
x2 = [xi_1+hw_ext,xi_1+hw_ext]
x3 = [xi_1-hw_ext,xi_1-hw_ext]

ax1.legend(loc=0)
for ax in [ax1,ax2]:
    ax.set_xlabel('$\\xi_3$ (eV)')
    ax.set_xlim([-1.,1.])


    ymin, ymax = ax.set_ylim()

    for x in [x1,x2,x3]:
        ax.plot(x,[ymin,ymax],'k--',label='__nolabel__')
    

ax1.set_ylabel('Re[$J$]')
ax2.set_ylabel('Im[$J$]')


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
