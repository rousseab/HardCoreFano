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

def generate_list_wedges(irreducible_wedge):

    list_wedges = []

    list_D6  = get_symmetries()
    for R in list_D6:
        list_k = []
        for k in irreducible_wedge.list_k:
            list_k.append(N.dot(R,k))
        list_k = N.array(list_k)

        # Make sure the Jacobian is always positive!
        triangles_indices = deepcopy(irreducible_wedge.triangles_indices)
        if N.linalg.det(R) < 0.:
            triangles_indices[:,1] = irreducible_wedge.triangles_indices[:,2]
            triangles_indices[:,2] = irreducible_wedge.triangles_indices[:,1]

        wedge = Wedge(list_k,triangles_indices)

        list_wedges.append(wedge)

    return list_wedges

def get_symmetries():
    theta = N.pi/3.
    # Identity
    E =  N.array([[  1.,  0.],
              [  0.,  1.]])


    # pi/3 rotation left
    C6_1  = N.array([[ N.cos(theta), N.sin(theta)],
             [-N.sin(theta), N.cos(theta)] ])

    # pi/3 rotation right
    C6_2  = N.array([[ N.cos(theta),-N.sin(theta)],
             [ N.sin(theta), N.cos(theta)] ])

    # 2pi/3 rotation left
    C3_1  = N.dot(C6_1,C6_1)

    # 2pi/3 rotation right
    C3_2  = N.dot(C6_2,C6_2)

    # pi rotation 
    C2    = N.dot(C6_1,C3_1)

    # vertical mirror plane, and rotations by 2pi/3
    # these mirror planes cross the face of the hexagon
    C2p_1  =  N.array([[ -1.,  0.],
               [  0.,  1.]])
    C2p_2  =  N.dot(C3_1,C2p_1)
    C2p_3  =  N.dot(C3_2,C2p_1)

    # mirror planes
    # these mirror planes cross the corners of the hexagon
    C2pp_1 =  N.dot(C6_1,C2p_1)
    C2pp_2 =  N.dot(C3_1,C2pp_1)
    C2pp_3 =  N.dot(C3_2,C2pp_1)


    list_D6 = N.array([ E, C6_1, C6_2, C3_1, C3_2, C2, C2p_1, C2p_2, C2p_3, C2pp_1, C2pp_2, C2pp_3])

    return list_D6 


def generate_irreducible_Wedge(nx_gridsize,fraction):
    """
    Generate the points inside a regular square
    """

    fy = fraction*N.sqrt(3.)/3*twopia # make the size just a tiny bit smaller, to avoid edge of 1BZ
    fx = fy/N.sqrt(3.)

    integers = N.arange(nx_gridsize+1)

    list_triangles = []

    ki = []
    kj = []

    triangles_indices  = []

    # create point array
    for j in integers:
        for i in integers[:j+1]:
            ki.append(i)
            kj.append(j)

    list_k = N.array([N.array(ki)*fx/nx_gridsize, N.array(kj)*fy/nx_gridsize]).transpose()

    # create connections
    for i in integers[:-1]:

        j  = i
        I1 = get_index(i,j+1)
        I2 = get_index(i,j)
        I3 = get_index(i+1,j+1)
        triangles_indices.append([I1,I2,I3])

        for j in integers[i+1:-1]:

            I1 = get_index(i,j+1)
            I2 = get_index(i,j)
            I3 = get_index(i+1,j+1)
            triangles_indices.append([I1,I2,I3])

            I1 = get_index(i+1,j)
            I2 = get_index(i+1,j+1)
            I3 = get_index(i,j)

            triangles_indices.append([I1,I2,I3])



    triangles_indices = N.array(triangles_indices )

    return list_k, triangles_indices



nmax_fine   = 8
nmax_coarse = 16


fig = plt.figure(figsize=(10,10))

ax  = fig.add_subplot(111)
ax.set_xlabel('$k_x$ ($2\pi/a$)')
ax.set_ylabel('$k_y$ ($2\pi/a$)')


fraction = 0.25
list_k_1, triangles_indices_1 = generate_irreducible_Wedge(nmax_fine,fraction)

fy = N.sqrt(3.)/3*twopia 
fx = fy/N.sqrt(3.)

dx = (1-fraction)*fx
dy = (1-fraction)*fy

translation = N.array([dx,dy])

list_k_2 = list_k_1+ translation 
triangles_indices_2 = triangles_indices_1 


list_k_fine = N.concatenate([list_k_1,list_k_2])
triangles_indices_fine = N.concatenate([triangles_indices_1, triangles_indices_2+len(list_k_1)])



list_k_3, triangles_indices_3 = generate_irreducible_Wedge(nmax_coarse,1.0)

# move the triangles around

dic_indices_coarse = {}
list_k_coarse = []

ic = -1

dx1 = fraction*fx
dy1 = fraction*fy

dx2 = (1-fraction)*fx
dy2 = (1-fraction)*fy

for i, k in enumerate(list_k_3):

    for1 = k[0] < dx1 and k[1] < dy1
    for2 = k[0] > dx2 and k[1] > dy2

    forbidden = for1 or for2

    if not forbidden:
        list_k_coarse.append(k)
        ic += 1
        dic_indices_coarse[i] = ic
    
list_k_coarse = N.array(list_k_coarse )

triangles_indices_coarse = []

for i1,i2,i3 in triangles_indices_3: 


    if dic_indices_coarse.has_key(i1) and dic_indices_coarse.has_key(i2) and dic_indices_coarse.has_key(i3):

        ic1 = dic_indices_coarse[i1]
        ic2 = dic_indices_coarse[i2]
        ic3 = dic_indices_coarse[i3]

        t = N.array([ic1,ic2,ic3])

        triangles_indices_coarse.append(t)

triangles_indices_coarse = N.array(triangles_indices_coarse )


list_k = N.concatenate([list_k_fine,list_k_coarse])
triangles_indices = N.concatenate([triangles_indices_fine, triangles_indices_coarse+len(list_k_fine)])


irreducible_wedge = Wedge(list_k,triangles_indices)


list_wedges = generate_list_wedges(irreducible_wedge)

for wedge in list_wedges:
    x = wedge.list_k[:,0]/twopia
    y = wedge.list_k[:,1]/twopia
    t = wedge.triangles_indices
    ax.triplot(x,y,triangles=t)


"""
x = list_k[:,0]/twopia
y = list_k[:,1]/twopia
t = triangles_indices
ax.triplot(x,y,triangles=t)
"""


"""
x = list_k_fine[:,0]/twopia
y = list_k_fine[:,1]/twopia
t = triangles_indices_fine
ax.triplot(x,y,triangles=t)


x = list_k_coarse[:,0]/twopia
y = list_k_coarse[:,1]/twopia
t = triangles_indices_coarse
ax.triplot(x,y,triangles=t)
"""




#x = list_k_1[:,0]/twopia
#y = list_k_1[:,1]/twopia
#t = triangles_indices_1
#ax.triplot(x,y,triangles=t)

#x = list_k_2[:,0]/twopia
#y = list_k_2[:,1]/twopia
#t = triangles_indices_2
#ax.triplot(x,y,triangles=t)



#x = list_k_3[:,0]/twopia
#y = list_k_3[:,1]/twopia
#t = triangles_indices_3
#ax.triplot(x,y,triangles=t)





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
