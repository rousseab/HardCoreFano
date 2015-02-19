import numpy as N
class Path:
    """
    This class generates and stores a kpath.
    """

    def __init__(self,acell,kpath_nsize,order=['M','$\Gamma$','K','M']):

        self.acell  = acell
        self.twopia = 2.*N.pi/acell

        self.G = N.array([0.,0.])
        self.M = N.array([N.sqrt(3.)/2.,0.5])*N.sqrt(3.)/3.*self.twopia 
        self.K = N.array([2./3.,0.])*self.twopia 

        self.special_points_dict = {'M':self.M,'$\Gamma$':self.G,'K':self.K}

        self.order = order

        self.generate_list_k(kpath_nsize)

    def generate_list_k(self,kpath_nsize):

        # Find the path length        
        path_length = 0.

        for L2,L1 in zip(self.order[1:],self.order[:-1]):

            P2  = self.special_points_dict[L2]
            P1  = self.special_points_dict[L1]
            L12 = N.linalg.norm(P2-P1)

            path_length += L12

        dx = path_length/kpath_nsize

        # Iterate on all labels
        self.list_x      = []
        self.list_k      = []
        self.list_xticks = []
        self.list_labels = []

        x  = 0.
        D  = 0.
        for L2,L1 in zip(self.order[1:],self.order[:-1]):

            P2  = self.special_points_dict[L2]
            P1  = self.special_points_dict[L1]
            L12 = N.linalg.norm(P2-P1)
            D  += L12
            dir_12 = (P2-P1)/L12

            k = P1
            self.list_x.append(x)
            self.list_k.append(k)
            self.list_xticks.append(x)
            self.list_labels.append(L1)

            while x < D-dx:
                x = x + dx
                k = k+ dir_12*dx

                self.list_x.append(x)
                self.list_k.append(k)

        
        self.list_xticks.append(x)
        self.list_labels.append(self.order[-1])
        self.list_x = N.array(self.list_x)
        self.list_k = N.array(self.list_k)
