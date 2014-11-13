#================================================================================
# Plot a nice tesselation of the 1BZ
#
#================================================================================
import common
reload(common)
from common import *

#----------------------------------------
# import modules
#----------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import  FontProperties

import matplotlib.cm as cm

mpl.rcParams['font.size'] = 20.
legendfonts = FontProperties(size=16)


def extract_unique_q_vectors(Q_grid):

    tol = 1e-12
    list_all_q      = []
    list_reduced_q  = []

    list_wedge_to_list = []
    for wedge in Q_grid.list_wedges:

        wedge_to_list = {}

        for ik, k in enumerate(wedge.list_k):

            list_all_q.append(k)
            new = True

            for jq, q in enumerate(list_reduced_q):

                test = N.linalg.norm(k-q) < tol
                if test:
                    wedge_to_list[ik] = jq
                    new = False
                    break

            if new:

                list_reduced_q.append(k)
                wedge_to_list[ik] = len(list_reduced_q)-1

        list_wedge_to_list.append(wedge_to_list)

    list_all_q      = N.array(list_all_q)
    list_reduced_q  = N.array(list_reduced_q)

    return list_all_q, list_reduced_q, list_wedge_to_list

nmax_coarse = 4
nmax_fine   = 8
n_blocks_coarse_to_fine = 1
include_Gamma = True
grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma )



list_all_q, list_reduced_q, list_wedge_to_list= extract_unique_q_vectors(grid)

# Test!
tol = 1e-12
for wedge, wtl in zip(grid.list_wedges,list_wedge_to_list):

    for ik, k in enumerate(wedge.list_k):
        q = list_reduced_q[wtl[ik]]
        error = N.linalg.norm(k-q)
        if error > tol:
            print 'wedge to list ERROR!'


fig = plt.figure(figsize=(10,10))

ax  = fig.add_subplot(111)
ax.set_xlabel('$k_x$ ($2\pi/a$)')
ax.set_ylabel('$k_y$ ($2\pi/a$)')


list_images = []

list_min = []
list_max = []

list_k = []




for i,wedge in enumerate(grid.list_wedges):

    x = wedge.list_k[:,0]/twopia
    y = wedge.list_k[:,1]/twopia
    t = wedge.triangles_indices

    F  = x**2+y**2
    fc = F[t].mean(axis=1)
    image = ax.tripcolor(x,y,triangles=t,facecolors=fc,cmap=plt.cm.spectral_r,  edgecolors='w')
    list_images.append(image)

    list_min.append(N.min(fc))
    list_max.append(N.max(fc))


vmin = N.min(list_min)
vmax = N.max(list_max)

norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
for image in list_images:
        image.set_norm(norm)

fig.colorbar(image,ax=ax)

for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
    line.set_markeredgewidth(3)


ax.set_xlim([-0.7,0.7])
ax.set_ylim([-0.7,0.7])
fig.gca().set_aspect('equal')
fig.suptitle('First Brillouin Zone Tesselation')


fig.subplots_adjust(    left    =       0.20,
                        bottom  =       0.20,
                        right   =       0.90,
                        top     =       0.90,
                        wspace  =       0.20,
                        hspace  =       0.20)

plt.show()
