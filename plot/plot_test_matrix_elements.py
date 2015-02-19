#================================================================================
# 
#       This script compares the computed matrix elements with the cone approximation
#       expectation for their values.
#
#
#================================================================================
import common
reload(common)
from common import *


import time
#----------------------------------------
# import modules
#----------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties

import matplotlib.cm as cm
import matplotlib.ticker as ticker

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

mpl.rcParams['font.size'] = 20.
legendfonts = FontProperties(size=16)

nmax_coarse = 8
nmax_fine   = 128
n_blocks_coarse_to_fine = 2
include_Gamma = False

clip_energy = 1.0

grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma,\
                                clip_grid=True,clip_energy=clip_energy)


c = N.cos(N.pi/3.)
s = N.sin(N.pi/3.)
list_K_points =twopia*2./3.*N.array([[ 1., 0.],  # K1
                                     [-c , s ],  # 
                                     [-c ,-s ],  # 
                                     [-1., 0.],  # K2
                                     [ c ,-s ],  # 
                                     [ c , s ]]) # 

K1_set    = [0,1,2]
K2_set    = [3,4,5]





path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_4.0/modules/'
list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)


#q_vector = N.random.random(2)/10.*twopia
q_vector = N.array([0.,0.])
list_hw, list_hw_renormalized, Eq = Renormalizor.get_frequency_and_polarization(q_vector)
nu = N.random.random_integers(0,5)
E_phonon_polarization = Eq[:,nu]
hw_nu_q               = list_hw[nu]


list_cmax_re = []
list_cmin_re = []
list_cmax_im = []
list_cmin_im = []

list_images_re = []
list_images_im = []


fig1 = plt.figure(figsize=(14,10))
fig2 = plt.figure(figsize=(14,10))

ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
list_ax = [ax1,ax2]
for ax in list_ax:
    ax.set_xlabel('$k_x$ ($2\pi/a$)')
    ax.set_ylabel('$k_y$ ($2\pi/a$)')


for i,wedge in enumerate(grid.list_wedges):

    representing_k = wedge.list_k[-1]

    d = N.sum((list_K_points-representing_k)**2,axis=1)
    imin = N.argmin(d)

    K_point = list_K_points[imin]    
    if imin in K1_set:
        eta = 1.
    elif imin in K2_set:
        eta =-1.

    list_cone_k  = wedge.list_k-K_point
    list_cone_nk = N.sqrt(N.sum(list_cone_k**2 ,axis=1))
    
    
    # The function fk        
    fk_cone = -N.sqrt(3.)/2.*acell*eta*(list_cone_k[:,0]+1j*eta*list_cone_k[:,1])
    fk        = function_fk(wedge.list_k)
    Fk = fk_cone - fk

    # The function Vk        
    Vkx, Vky   = function_Vk(wedge.list_k)

    Vkx_cone   = 1.5*1j*eta*N.ones_like(Vkx)
    Vky_cone   =-1.5*N.ones_like(Vky)

    #Fk = Vkx-Vkx_cone
    Fk = N.abs(Vkx-Vkx_cone)**2+N.abs(Vky-Vky_cone)**2


    # The Current
    M = MGridFunction(q_vector,E_phonon_polarization ,hw_nu_q, wedge)

    J_x_cone_11 =  hvF/Ha_to_eV*list_cone_k[:,0]/list_cone_nk 
    J_x_cone_22 = -hvF/Ha_to_eV*list_cone_k[:,0]/list_cone_nk 
    J_x_cone_12 = -1j*eta*hvF/Ha_to_eV*list_cone_k[:,1]/list_cone_nk 
    J_x_cone_21 =  1j*eta*hvF/Ha_to_eV*list_cone_k[:,1]/list_cone_nk 

    J_y_cone_11 =  hvF/Ha_to_eV*list_cone_k[:,1]/list_cone_nk 
    J_y_cone_22 = -hvF/Ha_to_eV*list_cone_k[:,1]/list_cone_nk 
    J_y_cone_12 =  1j*eta*hvF/Ha_to_eV*list_cone_k[:,0]/list_cone_nk 
    J_y_cone_21 = -1j*eta*hvF/Ha_to_eV*list_cone_k[:,0]/list_cone_nk 


    J_x_cone = N.array([[J_x_cone_11 ,J_x_cone_12],[J_x_cone_21,J_x_cone_22]])   
    J_y_cone = N.array([[J_y_cone_11 ,J_y_cone_12],[J_y_cone_21,J_y_cone_22]])   

    error_J_x = M.J_x-J_x_cone
    error_J_y = M.J_y-J_y_cone

    Fk        = error_J_y[1,1,:] 




    x = wedge.list_k[:,0]/twopia
    y = wedge.list_k[:,1]/twopia
    t = wedge.triangles_indices
    fc = (i+1)*N.arange(len(t))[::-1]
        
    image_re = ax1.tripcolor(x,y, N.real(Fk), shading='gouraud', edgecolors='k',cmap=plt.cm.spectral_r)
    image_im = ax2.tripcolor(x,y, N.imag(Fk), shading='gouraud', edgecolors='k',cmap=plt.cm.spectral_r)

    list_images_re.append(image_re)
    list_images_im.append(image_im)

    cmin_re, cmax_re = image_re.get_clim()
    cmin_im, cmax_im = image_im.get_clim()

    list_cmax_re.append(cmax_re)
    list_cmin_re.append(cmin_re)

    list_cmax_im.append(cmax_im)
    list_cmin_im.append(cmin_im)



cmin_re = N.min(list_cmin_re)
cmax_re = N.max(list_cmax_re)

cmin_im = N.min(list_cmin_im)
cmax_im = N.max(list_cmax_im)

for image in list_images_re:
    image.set_clim([cmin_re,cmax_re])

axc_re = fig1.colorbar(image,format=ticker.FuncFormatter(fmt))


for image in list_images_im:
    image.set_clim([cmin_im,cmax_im])

axc_im = fig2.colorbar(image,format=ticker.FuncFormatter(fmt))


#axc_re.set_label(r'$a_0^2/Ha$')
#axc_im.set_label(r'$a_0^2/Ha$')

# Plot Fermi surface

for ax in list_ax:
    for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
        line.set_markeredgewidth(3)


    ax.set_xlim([-0.7,0.7])
    ax.set_ylim([-0.7,0.7])
    fig1.gca().set_aspect('equal')
    fig2.gca().set_aspect('equal')


title = 'Singular $I$, q = (%3.2f, %3.2f) 2 $\pi/a$'%(q_vector[0]/twopia,q_vector[1]/twopia)

fig1.suptitle(title+'\n Real Part')
fig2.suptitle(title+'\n Imaginary Part')

for fig in [fig1,fig2]:
    fig.subplots_adjust(    left    =       0.10,
                            bottom  =       0.10,
                            right   =       0.90,
                            top     =       0.90,
                            wspace  =       0.20,
                            hspace  =       0.20)

plt.show()
