import numpy as N
class Path:
    """
    This class generates and stores a kpath.
    """

    def __init__(self,acell,kpath_nsize):

        self.acell  = acell
        self.twopia = 2.*N.pi/acell

        self.G = N.array([0.,0.])
        self.M = N.array([N.sqrt(3.)/2.,0.5])*N.sqrt(3.)/3.*self.twopia 
        self.K = N.array([2./3.,0.])*self.twopia 

        self.generate_list_k(kpath_nsize)

    def generate_list_k(self,kpath_nsize):

        LGK = N.linalg.norm(self.K-self.G)
        LKM = N.linalg.norm(self.M-self.K)
        LMG = N.linalg.norm(self.G-self.M)

        dir_GK = (self.K-self.G)/LGK
        dir_KM = (self.M-self.K)/LKM
        dir_MG = (self.G-self.M)/LMG


        path_length = LGK+LKM+LMG

        dx = path_length/kpath_nsize

        x  = 0.
        k  = 1.*self.G
        self.list_x = [x]
        self.list_k = [k]
        self.list_xticks = [x]
        self.list_labels = ['$\Gamma$']

        while x < LGK-dx:
            x = x + dx
            k = k+ dir_GK*dx

            self.list_x.append(x)
            self.list_k.append(k)

        x = 1.*LGK
        k = 1.*self.K
        self.list_x.append(x)
        self.list_k.append(k)

        self.list_xticks.append(x)
        self.list_labels.append('K')

        while x < LGK+LKM-dx:
            x = x+dx
            k = k+dir_KM*dx

            self.list_x.append(x)
            self.list_k.append(k)

        x = LGK+LKM
        k = 1.*self.M
        self.list_x.append(x)
        self.list_k.append(k)
        self.list_xticks.append(x)
        self.list_labels.append('M')


        while x < path_length-dx:
            x = x+dx
            k = k+dir_MG*dx
            self.list_x.append(x)
            self.list_k.append(k)

        x = 1.*path_length
        k = 1.*self.G
        self.list_x.append(x)
        self.list_k.append(k)
        self.list_xticks.append(x)
        self.list_labels.append('$\Gamma$')
        
        self.list_x = N.array(self.list_x)
        self.list_k = N.array(self.list_k)
