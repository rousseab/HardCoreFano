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

nmax = 4
grid = TesselationGrid(nmax)


fig = plt.figure(figsize=(10,10))

ax  = fig.add_subplot(111)
ax.set_xlabel('$k_x$ ($2\pi/a$)')
ax.set_ylabel('$k_y$ ($2\pi/a$)')


Integral = 0.
for i,wedge in enumerate(grid.list_wedges):

	ones = N.ones(len(wedge.list_k))
	list_Fk = ones[:,N.newaxis]

	Integral += AreaIntegrator(wedge,list_Fk)[0]

	x = wedge.list_k[:,0]/twopia
	y = wedge.list_k[:,1]/twopia
	t = wedge.triangles_indices
	fc = (i+1)*N.arange(len(t))[::-1]
        

	ax.triplot(x,y,triangles=t)


FBZ_area = (2.*N.pi)**2/Area
error = FBZ_area - Integral

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
