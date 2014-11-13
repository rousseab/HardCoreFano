#================================================================================
#
#                  Module Grids
#           ===========================
#
#   The various classes implemented in this module are used to compute
#   sums on the first Brillouin zone by creating grids and tesselations.
#================================================================================

from module_Constants import *
import sys

# utility function
def get_index(i,j):
    I = j*(j+1)/2+i
    return I

class Wedge:
    """
    This is a container class which stores the 
    k points and their connection topology to form
    an irreducible wedge of the first Brillouin zone.
    """

    def __init__(self,list_k,triangles_indices):
        self.list_k            = list_k
        self.triangles_indices = triangles_indices

        # create cross product for every triangle

        self.cross_products = []
        self.Jacobians      = []

        for t in self.triangles_indices:
            k1 = self.list_k[t[0]]
            k2 = self.list_k[t[1]]
            k3 = self.list_k[t[2]]

            self.cross_products.append( N.linalg.norm(N.cross(k2-k1,k3-k1)) )

            J = complex(1.,0.)*N.array([ [k2[0]-k1[0],k3[0]-k1[0]] ,\
                                         [k2[1]-k1[1],k3[1]-k1[1]] ])
        
    
            self.Jacobians.append(J)


        self.cross_products = N.array(self.cross_products)
        self.Jacobians      = N.array(self.Jacobians)


class TesselationGrid:
    """
    This class generates and stores a k-grid which is consistent
    with the hexagonal symmetry of graphene.
    """

    def __init__(self,nx_gridsize):

        self.get_symmetries()

        delta = 0.001
        fy = (N.sqrt(3.)/3.-delta)*twopia # make the size just a tiny bit smaller, to avoid edge of 1BZ
        fx = fy/N.sqrt(3.)

        self.irreducible_wedge = self.generate_irreducible_Wedge(nx_gridsize,fx,fy)
        self.generate_list_wedges()

    def generate_list_wedges(self):

        self.list_wedges = []

        self.irreducible_wedge.triangles_indices
        for R in self.list_D6:
            list_k = []
            for k in self.irreducible_wedge.list_k:
                list_k.append(N.dot(R,k))

            list_k = N.array(list_k)

            # Make sure the Jacobian is always positive!
            triangles_indices = deepcopy(self.irreducible_wedge.triangles_indices)
            if N.linalg.det(R) < 0.:
                triangles_indices[:,1] = self.irreducible_wedge.triangles_indices[:,2]
                triangles_indices[:,2] = self.irreducible_wedge.triangles_indices[:,1]

            wedge = Wedge(list_k,triangles_indices)

            self.list_wedges.append(wedge)



    def generate_irreducible_Wedge(self,nx_gridsize,fx,fy):
        """
        Generate the points inside a regular square
        """

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
        irreducible_wedge = Wedge(list_k=list_k,triangles_indices=triangles_indices)

        return irreducible_wedge 


    def get_symmetries(self):
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


        self.list_D6 = N.array([ E, C6_1, C6_2, C3_1, C3_2, C2, C2p_1, C2p_2, C2p_3, C2pp_1, C2pp_2, C2pp_3])



class TesselationDoubleGrid(TesselationGrid):
    """
    This class generates and stores a k-grid which is consistent
    with the hexagonal symmetry of graphene.

    A fine grid around Gamma and K, and a coarse grid in the rest of the 
    zone are generated.
    """

    def __init__(self,nx_gridsize_coarse, nx_gridsize_fine, n_blocks_coarse_to_fine, include_Gamma):


        # Test various parameters to make sure they are sound

        if nx_gridsize_fine%nx_gridsize_coarse != 0:
            print 'ERROR! The fine and coarse grids must be compatible! make sure nx_gridsize_fine = (integer) x nx_gridsize_coarse'
            sys.exit()              

        if n_blocks_coarse_to_fine < 1 or n_blocks_coarse_to_fine > nx_gridsize_coarse:
            print 'ERROR! Pick a reasonable amount of blocks for the fine grid!'
            sys.exit()              

        self.fraction = (1.*n_blocks_coarse_to_fine)/(1.*nx_gridsize_coarse)

        self.get_symmetries()

        delta = 0.0001
        self.fy = (N.sqrt(3.)/3.-delta)*twopia # make the size just a tiny bit smaller, to avoid edge of 1BZ
        self.fx = self.fy/N.sqrt(3.)

        self.include_Gamma = include_Gamma

        self.nx_gridsize_coarse = nx_gridsize_coarse  
        self.nx_gridsize_fine   = self.fraction*nx_gridsize_fine

        self.irreducible_wedge = self.generate_irreducible_Wedge_double_grid()


        self.generate_list_wedges()


    def generate_irreducible_Wedge_double_grid(self):


        # First, generate a dense grid in the bottom corner of the wedge
        small_fx = self.fraction*self.fx
        small_fy = self.fraction*self.fy

        dummy_wedge = self.generate_irreducible_Wedge(self.nx_gridsize_fine,small_fx,small_fy)

        list_k_1 = dummy_wedge.list_k
        triangles_indices_1 = dummy_wedge.triangles_indices

        # Translate the small triangle to the top right corner
        dx = (1.-self.fraction)*self.fx
        dy = (1.-self.fraction)*self.fy

        translation = N.array([dx,dy])

        list_k_2 = list_k_1+ translation 
        triangles_indices_2 = deepcopy(triangles_indices_1)


        if self.include_Gamma: 
            # Create the list of k points for the fine sub-grid
            list_k_fine = N.concatenate([list_k_1,list_k_2])
            triangles_indices_fine = N.concatenate([triangles_indices_1, triangles_indices_2+len(list_k_1)])
        else:
            list_k_fine = list_k_2
            triangles_indices_fine = triangles_indices_2


        # Create the coarse grid
        dummy_wedge = self.generate_irreducible_Wedge(self.nx_gridsize_coarse,self.fx,self.fy)
        list_k_3 = dummy_wedge.list_k
        triangles_indices_3 = dummy_wedge.triangles_indices


        # move the triangles around
        dic_indices_coarse = {}
        list_k_coarse = []

        ic = -1

        dx1 = self.fraction*self.fx
        dy1 = self.fraction*self.fy

        dx2 = (1.-self.fraction)*self.fx
        dy2 = (1.-self.fraction)*self.fy

        tol = 1e-8

        for i, k in enumerate(list_k_3):


            if self.include_Gamma: 
                for1 = k[0] < dx1-tol and k[1] < dy1-tol
                for2 = k[0] > dx2+tol and k[1] > dy2+tol
                forbidden = for1 or for2
            else: 
                forbidden = k[0] > dx2+tol and k[1] > dy2+tol

            #for1 = k[0] < dx1 and k[1] < dy1
            #for2 = k[0] > dx2 and k[1] > dy2


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

        triangles_indices = N.concatenate([triangles_indices_fine, triangles_indices_coarse+len(list_k_fine)]).astype(N.int)


        irreducible_wedge = Wedge(list_k,triangles_indices)

        return irreducible_wedge 

