#================================================================================
# Plot a nice tesselation of the 1BZ
#
#================================================================================
import sys
module_directory = '/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_4.0/modules'



sys.path.append(module_directory)

from module_Driver import *


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties
import matplotlib.cm as cm
mpl.rcParams['font.size'] = 32.
legendfonts = FontProperties(size=24)

from Scientific.IO.NetCDF import NetCDFFile as Dataset

mpl.rcParams['font.size'] = 40.

import matplotlib.tri as tri
import warnings
warnings.filterwarnings("ignore")

#nq_max = 64
nq_max = 32
nk_max = 64
mu     = -0.400
amma  = 0.050
#filename = './HCF_nq_max=%i_nk_max=%i_mu=%i_meV_Gamma=%i_meV.nc'%(nq_max,nk_max,1000*mu,1000*Gamma)

filename = './HCF_nq_max=%i_nk_max=%i_mu=%i_meV_Gamma=%i_meV.nc'%(nq_max,nk_max,1000*mu,1000*Gamma)
print filename

SCALE = 1e4

ihw       = N.int(150)
#ihw       = N.int(0)


ncfile   = Dataset(filename,'r')
list_hw  = ncfile.variables['list_hw'][:]

hw       = list_hw[ihw]

k_mu     = N.abs(mu)/hvF
k_w      = N.abs(hw)/hvF



list_q_plus    = ncfile.variables['list_q_plus'][:]
list_q_minus   = ncfile.variables['list_q_minus'][:]
list_q         = N.concatenate([list_q_plus,list_q_minus])


list_nu        = ncfile.variables['phonon_nu'][:]


list_S_xx = N.zeros_like(list_q[:,0])
list_S_yy = N.zeros_like(list_q[:,0])


list_inu = [0,1,2,3]


u_const = u_a0
v_const = 0.


for inuu in list_inu:

    inu = N.int(inuu)
    nu  = list_nu[inu]

    print 'nu = %i'%nu
    phonon_frequencies = ncfile.variables[ 'phonon_frequencies'][inuu,:]
    print phonon_frequencies

    H_plus  = u_const*(complex(1.,0.)*ncfile.variables[ 'Re_H_plus'][inu,:,0,:,ihw]  +\
                       complex(0.,1.)*ncfile.variables[ 'Im_H_plus'][inu,:,0,:,ihw]) +\
              v_const*(complex(1.,0.)*ncfile.variables[ 'Re_H_plus'][inu,:,1,:,ihw]  +\
                       complex(0.,1.)*ncfile.variables[ 'Im_H_plus'][inu,:,1,:,ihw])

    H_minus = u_const*(complex(1.,0.)*ncfile.variables[ 'Re_H_minus'][inu,:,0,:,ihw]  +\
                       complex(0.,1.)*ncfile.variables[ 'Im_H_minus'][inu,:,0,:,ihw]) +\
              v_const*(complex(1.,0.)*ncfile.variables[ 'Re_H_minus'][inu,:,1,:,ihw]  +\
                       complex(0.,1.)*ncfile.variables[ 'Im_H_minus'][inu,:,1,:,ihw])


    hqhq_xx   = N.abs(H_plus[:,0]*H_minus[:,0])
    hqhq_yy   = N.abs(H_plus[:,1]*H_minus[:,1])

    list_S_xx += N.concatenate([hqhq_xx, hqhq_xx])
    list_S_yy += N.concatenate([hqhq_yy, hqhq_yy])


ncfile.close()



list_titles = [r'$|\sum_{\nu} H^{\nu x}({{\bf q},\hbar\omega})\times H^{\nu x}({{-\bf q},\hbar\omega})$|, $\hbar\omega$ = %i meV'%(1000*hw),
               r'$|\sum_{\nu} H^{\nu y}({{\bf q},\hbar\omega})\times H^{\nu y}({{-\bf q},\hbar\omega})$|, $\hbar\omega$ = %i meV'%(1000*hw)]

list_names = ['HxHx.pdf','HyHy.pdf']
list_alpha = ['x','y']

c = N.cos(N.pi/3.)
s = N.sin(N.pi/3.)
list_centers = N.array([ [0.,0.], [-2./3.,0.], [ 2./3.,0.], 2./3.*N.array([c,s]), 
                       2./3.*N.array([c,-s]), 2./3.*N.array([-c,s]), 2./3.*N.array([-c,-s])]) 


list_datafilenames = ['Sxx.dat','Syy.dat']

for datafilename, list_S, title, name, alpha in zip(list_datafilenames, [list_S_xx,list_S_yy],list_titles, list_names,list_alpha):


    file = open(datafilename,'w')

    print >> file, '#================================================================================'
    print >> file, '#'
    print >> file, '#       S function, within the first Brillouin zone'
    print >> file, '#'
    print >> file, '#    qx (2pi/a)       qy (2pi/a)      S (10^4 e hbar/m a0)'
    print >> file, '#================================================================================'


    x = list_q[:,0]/twopia
    y = list_q[:,1]/twopia


    for xx, yy, SS in zip(x,y,SCALE*list_S):
        print >> file, '%16.12f  %16.12f  %16.12f'%(xx,yy,SS)

    file.close()


    fig = plt.figure(figsize=(15,12))

    ax  = fig.add_subplot(111)
    ax.set_xlabel('$q_x$ ($2\pi/a$)')
    ax.set_ylabel('$q_y$ ($2\pi/a$)')


    mymap = plt.get_cmap("spectral")

    colors = mymap(SCALE*list_S)

    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triang = tri.Triangulation(x, y) 

    image = ax.tripcolor(triang, SCALE*list_S, shading='gouraud', edgecolors='k',cmap=plt.cm.spectral_r)
    #image = ax.tripcolor(triang, list_S, shading='flat', edgecolors='k',cmap=plt.cm.spectral)
    axc = plt.colorbar(image)

    #axc.set_label(r'$S^{%s}({\bf q},\omega)$'%alpha)
    axc.set_label(r'$\times 10^{%+i}$ $(e \hbar /m a_0)^2$'%-4)

    #ax.scatter(x,y,s=10)
    rad_max = (2*k_mu+k_w)/twopia
    rad_mid = (2*k_mu-k_w)/twopia
    rad_min = k_w/twopia
    th = N.arange(0.,2.*N.pi,0.01)

    nb = 0
    for center in list_centers:

        x = center[0]+rad_max*N.cos(th)
        y = center[1]+rad_max*N.sin(th)

        ax.plot(x, y,'k-',lw=2)


        nb += 1
        circle_file = open('circle_%i.dat'%nb,'w')
        print >> circle_file,"#================================"
        print >> circle_file,"# Coordinates for Circle  %i"%nb
        print >> circle_file,"#================================"
        print >> circle_file,'#    qx (2pi/a)       qy (2pi/a)'
        print >> circle_file,"#================================"
        for xcart, ycart in zip(x,y):
            print >> circle_file, "%12.8f     %12.8f"%(xcart,ycart)

        circle_file.close() 


        #x = center[0]+rad_mid*N.cos(th)
        #y = center[1]+rad_mid*N.sin(th)
        #ax.plot(x, y,'k-',lw=2)

        x = center[0]+rad_min*N.cos(th)
        y = center[1]+rad_min*N.sin(th)
        ax.plot(x, y,'k-',lw=2)

        nb += 1
        circle_file = open('circle_%i.dat'%nb,'w')
        print >> circle_file,"#================================"
        print >> circle_file,"# Coordinates for Circle  %i"%nb
        print >> circle_file,"#================================"
        print >> circle_file,'#    qx (2pi/a)       qy (2pi/a)'
        print >> circle_file,"#================================"
        for xcart, ycart in zip(x,y):
            print >> circle_file, "%12.8f     %12.8f"%(xcart,ycart)

        circle_file.close() 




        for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
            line.set_markeredgewidth(3)


        ax.set_xlim([-0.7,0.7])
        ax.set_ylim([-0.7,0.7])
        fig.gca().set_aspect('equal')

        #fig.suptitle(title)


        fig.subplots_adjust(    left    =       0.16,
                                bottom  =       0.16,
                                right   =       0.96,
                                top     =       0.98,
                                wspace  =       0.20,
                                hspace  =       0.20)


    plt.savefig(name)
