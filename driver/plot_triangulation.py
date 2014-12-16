#================================================================================
# Plot a nice tesselation of the 1BZ
#
#================================================================================
import sys
module_directory = '/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_5.0/modules'

sys.path.append(module_directory)

from module_Driver import *


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties
import matplotlib.cm as cm
mpl.rcParams['font.size'] = 16.
legendfonts = FontProperties(size=16)

import matplotlib.ticker as ticker



from Scientific.IO.NetCDF import NetCDFFile as Dataset

mpl.rcParams['font.size'] = 16.

import matplotlib.tri as tri
import warnings
warnings.filterwarnings("ignore")



def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

mu     = -0.400
Gamma  = 0.050

filename = 'HCF_nq=8_32_3_nk=8_32_3_mu=-400_meV_Gamma=50_meV_delta=50_meV.nc'
print filename

ihw       = N.int(150)


ncfile   = Dataset(filename,'r')
list_hw  = ncfile.variables['list_hw'][:]

hw       = list_hw[ihw]
print 'hw = %8.4f'%hw

k_mu     = N.abs(mu)/hvF
k_w      = N.abs(hw)/hvF



list_q      = ncfile.variables['list_q'][:]
list_nu     = ncfile.variables['phonon_nu'][:]

list_S      = N.zeros_like(list_q[:,0])
list_S_re   = N.zeros_like(list_q[:,0])
list_S_im   = N.zeros_like(list_q[:,0])


list_inu = [0,1,2,3]


u_const = u_a0

for inuu in list_inu:

    inu = N.int(inuu)
    nu  = list_nu[inu]

    print 'nu = %i'%nu
    phonon_frequencies = ncfile.variables[ 'phonon_frequencies'][inuu,:]

    # hqhq dimensions: nu, number_of_q_points,  number_of_frequencies
    re_hqhq = ncfile.variables['Re_HH'][inu,:,ihw]
    im_hqhq = ncfile.variables['Im_HH'][inu,:,ihw]

    # Sum over u parameters
    hqhq    = u_const**2*(complex(1.,0.)*re_hqhq[:] + complex(0.,1.)*im_hqhq[:])

    list_S    += N.abs(hqhq)
    list_S_re += N.real(hqhq)
    list_S_im += N.imag(hqhq)

ncfile.close()



title = r'$|\sum_{\nu\alpha} H^{\nu\alpha }({{\bf q},\hbar\omega})\times H^{\nu \alpha}({{-\bf q},\hbar\omega})$|, $\hbar\omega$ = %i meV'%(1000*hw)

SS = [list_S,list_S_re,list_S_im]
list_names = ['HH_abs.pdf','HH_real.pdf','HH_im.pdf']

c = N.cos(N.pi/3.)
s = N.sin(N.pi/3.)
list_centers = N.array([ [0.,0.], [-2./3.,0.], [ 2./3.,0.], 2./3.*N.array([c,s]), 
                       2./3.*N.array([c,-s]), 2./3.*N.array([-c,s]), 2./3.*N.array([-c,-s])]) 

list_datafilenames = ['S.dat','S_re.dat','S_im.dat']

for datafilename, list_S, name in zip(list_datafilenames, SS, list_names):


    S_data = list_S

    file = open(datafilename,'w')

    print >> file, '#================================================================================'
    print >> file, '#'
    print >> file, '#       S function, within the first Brillouin zone'
    print >> file, '#'
    print >> file, '#    qx (2pi/a)       qy (2pi/a)      S (e hbar/m a0)'
    print >> file, '#================================================================================'


    x = list_q[:,0]/twopia
    y = list_q[:,1]/twopia


    for xx, yy, SS in zip(x,y,S_data):
        print >> file, '%16.12f  %16.12f  %22.16f'%(xx,yy,SS)

    file.close()



    fig = plt.figure(figsize=(15,12))

    ax  = fig.add_subplot(111)
    ax.set_xlabel('$q_x$ ($2\pi/a$)')
    ax.set_ylabel('$q_y$ ($2\pi/a$)')


    mymap = plt.get_cmap("spectral_r")


    colors = mymap(S_data)

    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triang = tri.Triangulation(x, y)

    image = ax.tripcolor(triang, S_data, shading='gouraud', edgecolors='k',cmap=plt.cm.spectral_r)
    #image = ax.tripcolor(triang, list_S, shading='flat', edgecolors='k',cmap=plt.cm.spectral)
    axc = plt.colorbar(image,format=ticker.FuncFormatter(fmt))


    #axc.set_label(r'$S^{%s}({\bf q},\omega)$'%alpha)
    axc.set_label(r'$(e \hbar /m a_0)^2$')

    #ax.scatter(x,y,s=10)
    rad_max = (2*k_mu+k_w)/twopia
    rad_mid = (2*k_mu-k_w)/twopia
    rad_min = k_w/twopia
    th = N.arange(0.,2.*N.pi,0.01)

    for center in list_centers:

        x = center[0]+rad_max*N.cos(th)
        y = center[1]+rad_max*N.sin(th)
        ax.plot(x, y,'k--',lw=1)

        x = center[0]+rad_mid*N.cos(th)
        y = center[1]+rad_mid*N.sin(th)
        ax.plot(x, y,'k--',lw=1)



        x = center[0]+rad_min*N.cos(th)
        y = center[1]+rad_min*N.sin(th)
        ax.plot(x, y,'k--',lw=1)



    for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
        line.set_markeredgewidth(3)


    ax.set_xlim([-0.7,0.7])
    ax.set_ylim([-0.7,0.7])
    fig.gca().set_aspect('equal')



    fig.subplots_adjust(    left    =       0.03,
                            bottom  =       0.16,
                            right   =       0.96,
                            top     =       0.98,
                            wspace  =       0.20,
                            hspace  =       0.20)


    plt.savefig(name)
