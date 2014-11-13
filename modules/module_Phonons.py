#================================================================================
#
#            Module Phonons 
#            ===========================
#
#    This module implements the force constants to extract the phonon
#    dispersions and polarizations.
#
#    This module implements a more brute force, but simpler, scheme.
#
#================================================================================
from module_Constants import *
from module_D6h import *


def get_R_matrix(q,list_FC,list_R):

    list_qR = N.dot(list_R,q)

    AmB     = N.exp(-1j*N.dot(tau_A-tau_B,q))
    BmA     = N.exp(-1j*N.dot(tau_B-tau_A,q))

    phases  = N.exp(-1j*list_qR)

    factor        = complex(1.,0.)*N.ones([6,6])
    factor[:3,3:] = factor[:3,3:]*AmB
    factor[3:,:3] = factor[3:,:3]*BmA


    nR = len(list_R)    
    list_phases = N.resize(phases,(6,6,nR)).transpose(2,0,1)*factor

    R_matrix = N.sum(list_phases*list_FC,axis=0)*Ha_to_eV

    return R_matrix

def get_L_matrix(q,list_FC,list_R):

    list_qR = N.dot(list_R,q)

    phases  = N.exp(-1j*list_qR)

    nR = len(list_R)    
    list_phases = N.resize(phases,(6,6,nR)).transpose(2,0,1)

    L_matrix = N.sum(list_phases*list_FC,axis=0)*Ha_to_eV

    return L_matrix

def get_frequencies_and_polarizations(q,list_FC,list_R):

    list_qR = N.dot(list_R,q)

    AmB     = N.exp(-1j*N.dot(tau_A-tau_B,q))
    BmA     = N.exp(-1j*N.dot(tau_B-tau_A,q))

    phases  = N.exp(-1j*list_qR)

    factor        = complex(1.,0.)*N.ones([6,6])
    factor[:3,3:] = factor[:3,3:]*AmB
    factor[3:,:3] = factor[3:,:3]*BmA


    nR = len(list_R)    
    list_phases = N.resize(phases,(6,6,nR)).transpose(2,0,1)*factor

    R_matrix = N.sum(list_phases*list_FC,axis=0)*Ha_to_eV

    eigs, Bq = N.linalg.eigh(R_matrix)

    Eq = 1.*Bq

    Eq[:3,:] = Eq[:3,:]*N.exp(1j*N.dot(tau_A,q))
    Eq[3:,:] = Eq[3:,:]*N.exp(1j*N.dot(tau_B,q))

    list_hw = N.sign(eigs)*N.sqrt(N.abs(eigs))

    return list_hw, Eq

class CompteMauriRenormalize:
    """
    This class will compute the renormalized (and unrenormalized) 
    phonon spectrum if q is within the cutoff radius.

    We make this a class so that certain data structures can be 
    generated on initialization, once and for all.
    """
    def __init__(self,mauri_filepath,list_FC,list_R):

        self.list_FC = list_FC
        self.list_R  = list_R

        # Prepare data which will be used to interpolate phonon branch
        data = N.loadtxt(mauri_filepath+'mauri_GW_K.txt')

        x  = data[:,0]
        xK = norm_K_point/(twopia)

        I1  = N.where(x <= xK)[0] # along the Gamma-K line
        I2  = N.where(x >= xK)[0] # along the K-M line

        x1   = xK-data[I1,0]
        hwM1 = data[I1,1]/eV_to_cm1 # in eV
        I    = N.argsort(x1)

        self.x1   = x1[I]
        self.hwM1 = hwM1[I]

        x2   = data[I2,0]-xK
        hwM2 = data[I2,1]/eV_to_cm1 # in eV
        I    = N.argsort(x2)

        self.x2   = x2[I]
        self.hwM2 = hwM2[I]

        self.get_fit_polynomials()

        # This cutoff is the approximate range of the smallest dataset
        self.Q_cutoff  = 4.5/30.*twopia

        # First valley
        K1 = norm_K_point*N.array([1.,0.])
        K2 = norm_K_point*N.array([-0.5, N.sqrt(3.)/2])
        K3 = norm_K_point*N.array([-0.5,-N.sqrt(3.)/2])

        self.list_Kv1 = N.array([K1,K2,K3])
        # second valley
        K4 = -K1
        K5 = -K2
        K6 = -K3
        list_hw_renormalized = []
        self.list_Kv2 = N.array([K4,K5,K6])


        self.prepareDisentanglement()



    def prepareDisentanglement(self):
        """
        Prepare some arrays necessary to perform the k.p perturbation
        theory, which will disentangle phonon branches in the vicinity of K.
        """

        self.list_P = N.concatenate([self.list_Kv1,self.list_Kv2])

        self.list_L_matrix = []
        self.list_E        = []
        self.list_eps      = []


        for P in self.list_P:

            L_matrix = get_L_matrix(P,self.list_FC,self.list_R)
            self.list_L_matrix.append(L_matrix)

            eigs, E  = N.linalg.eigh(L_matrix)

            # Keep highest branch only, which we know to be the A1' mode at K.
            self.list_E.append(E[:,5])
            self.list_eps.append(eigs[5])

        return 


    def get_fit_polynomials(self):
        """
        Let's fit Mauri's data with a polynomial (as opposed to a spline).
        This way, small errors due to manually extracting points from a digitized image
        will be washed away. Also, splines yield exactly zero the instant the argument
        exits the range: this leads to large kinks when I accidentally leave the range...
        """

        # Playing around suggests this works well
        deg   = 8

        self.pfit1 = N.polyfit(self.x1,self.hwM1,deg=deg)
        self.pfit2 = N.polyfit(self.x2,self.hwM2,deg=deg)

        return

    def get_r1(self,theta):
    
        pi_on_3 = N.pi/3.
        th = N.mod(theta,2*pi_on_3)/pi_on_3

        if th <= 1.:
            r1 = th
        else:
            r1 = 2.-th

        return r1

    def get_angle(self,vector):

        tol = 1e-12
        x = vector[0]
        y = vector[1]


        # Find the angle
        if N.abs(x) < tol and N.abs(y) < tol:
            theta = 0.
        elif N.abs(x) < tol:

            fac   = (1.-N.sign(y))/2.
            theta = N.pi/2. + fac*N.pi
        else:
            fac = 0.5*(1.-N.sign(x))

            theta = N.mod( fac*N.pi + N.arctan(y/x), 2.*N.pi )

        return theta

    def kdotp_perturbation(self,iK,E,L):

        L0   = self.list_L_matrix[iK]
        e0   = self.list_E[iK]
        eps0 = self.list_eps[iK]

        ea = E[:,5]
        eb = E[:,4]

        #----------------------------------
        # Implement linear response for E
        #----------------------------------

        # the "perturbation"
        dL = L-L0

        # the change in eigenvalue
        deps = N.dot( N.conjugate(e0),N.dot(dL,e0))

        # the right-hand-size of the eigenvalue equation, to linear order
        rhs = deps*e0-N.dot(dL,e0)

        # the solution dE
        M  = L0-N.identity(6)*eps0+0.1*e0[:,N.newaxis]*N.conjugate(e0)

        de = N.dot(N.linalg.inv(M),rhs)

        # approximate solution for top branch
        e_approx = e0+de
        e_approx = e_approx/N.linalg.norm(e_approx)


        # test which actual branch is closest to approx
        testa = 1.-N.abs(N.dot(N.conj(e_approx),ea))**2
        testb = 1.-N.abs(N.dot(N.conj(e_approx),eb))**2

        if testa < testb:
            renormalization_index = 5
        else:
            renormalization_index = 4


        return renormalization_index 

    def get_frequency_and_polarization(self,q):
        """
        Compute the frequencies and polarization vectors
        for this q, for the force constant model and for the 
        Mauri renormalization.
        """

        # First, the eigenvalues and eigenvectors for the 
        # force constant model
        list_qR = N.dot(self.list_R,q)

        phases  = N.exp(-1j*list_qR)

        nR = len(self.list_R)    

        list_phases = N.resize(phases,(6,6,nR)).transpose(2,0,1)

        L_matrix = N.sum(list_phases*self.list_FC,axis=0)*Ha_to_eV

        eigs, Eq = N.linalg.eigh(L_matrix)

        list_hw  = N.sign(eigs)*N.sqrt(N.abs(eigs))

        list_hw_renormalized  = 1.*list_hw  

        # Find if the spectrum must be renormalized

        # Find if q is within a valley
        renormalize = False

        iK = -1

        for K in self.list_Kv1:
            iK += 1
            v = q-K
            if N.linalg.norm(v) < self.Q_cutoff:
                renormalize = True
                Q = v
                break

        if not renormalize:
            for K in self.list_Kv2:
                iK += 1
                v = q-K
                if N.linalg.norm(v) < self.Q_cutoff:
                    renormalize = True
                    Q =-v
                    break


        if renormalize:

            # We are in a renormalization valley!

            # First, identify which frequency must, in fact, be renormalized
            renormalization_index = self.kdotp_perturbation(iK,Eq,L_matrix)

            #print iK, renormalization_index 

            hw0 = list_hw[renormalization_index]

            # second, find the renormalized frequency
            x = N.linalg.norm(Q)/twopia

            # Find fits along the G-K and K-M lines
            hw_fit1 = N.polyval(self.pfit1,x)
            hw_fit2 = N.polyval(self.pfit2,x)


            theta = self.get_angle(Q)
            r1    = self.get_r1(theta)
            r2    = 1.-r1

            hw_fit = r1*hw_fit1 + r2*hw_fit2

            d = N.linalg.norm(Q)/self.Q_cutoff

            hw_renormalized = (1.-d**4)*hw_fit+d**4*hw0

            list_hw_renormalized[renormalization_index] = hw_renormalized 

            if renormalization_index == 4:
                # inverse the order of the two top branches, to avoid nasty crossings in plots
                x = list_hw[4]
                y = list_hw[5]

                list_hw[5] = 1.*x
                list_hw[4] = 1.*y

                x = list_hw_renormalized[4]
                y = list_hw_renormalized[5]

                list_hw_renormalized[5] = 1.*x
                list_hw_renormalized[4] = 1.*y

                ex = Eq[:,4]
                ey = Eq[:,5]

                Eq[:,5] = 1.*ex
                Eq[:,4] = 1.*ey


        return list_hw, list_hw_renormalized, Eq



