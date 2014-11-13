#================================================================================
#
#           Module M
#           ===============
#
#       This module implements an object which represents the matrix M,
#       namely the quantity which does not depend on frequency in the 
#       formalism.
#
#================================================================================

from module_Constants import *
from module_MatrixList import *
from module_Grids import *
from module_Functions import *


class MGridFunction:
    """
    Class which builds and contains the M matrix, which essentially represents 
    all the vertices in the diagram, and is frequency independent.

    The computations are separated in a component which depends on the u and v 
    parameters, related to the electron-phonon coupling.
    """

    def __init__(self,q_vector,E_phonon_polarization,hw_nu_q, wedge):

        self.q_vector = q_vector
        self.E_ph     = E_phonon_polarization
        self.hw_nu_q  = hw_nu_q

        self.build_indices()

        self.nk   = len(wedge.list_k)
        self.dim  = len(self.index_dictionary)

        # M will be in units of [charge] x [velocity] x [energy],
        # where all terms are in fundamental units

        # the second dimension is for u and v
        self.M    = complex(0.,0.)*N.zeros([self.dim,2,self.nk])

        self.epsilon_k = complex(0.,0.)*N.zeros([2,self.nk])
        self.epsilon_kq= complex(0.,0.)*N.zeros([2,self.nk])

        self.build_M(wedge)


    def build_indices(self):

        self.index_dictionary = {}
            
        index = -1
        for alpha in ['x','y']:
            for L in ['A','B']:
                for n1 in [0,1]:
                    for n2 in [0,1]:
                        for n3 in [0,1]:
                            index += 1 
                            key = (alpha,L,n1,n2,n3)
                            self.index_dictionary[key] = index



    def build_M(self,wedge):
        """
        All quantities are built here. It will be *crucial* to make
        this code as transparent as possible, to avoid making mistakes!
        """

        list_k    = wedge.list_k
        list_kq   = list_k + self.q_vector

        fk        = function_fk(list_k)
        abs_fk    = N.abs(fk)
        Ak        = fk/abs_fk
        Ak_star   = N.conjugate(Ak)
        eps_k     = function_epsilon_k(list_k)

        self.epsilon_k[0,:]  =-eps_k
        self.epsilon_k[1,:]  = eps_k

        fkq       = function_fk(list_kq)
        abs_fkq   = N.abs(fkq)
        Akq       = fkq/abs_fkq
        Akq_star  = N.conjugate(Akq)
        eps_kq    = function_epsilon_k(list_kq)

        self.epsilon_kq[0,:] =-eps_kq
        self.epsilon_kq[1,:] = eps_kq


        Vkx, Vky   = function_Vk(list_k)
        Vkqx, Vkqy = function_Vk(list_kq)

        Vqx, Vqy   = function_Vk(list_kq-list_k)

        # Build various Matrices 
        # Usage: MatrixList(m_11,m_12,m_21,m_22)

        ones  = complex(1.,0.)*N.ones(self.nk)
        zeros = complex(1.,0.)*N.zeros(self.nk)

        # Projectors on site A and site B
        PA   = MatrixList(    ones ,   zeros,\
                              zeros,   zeros)

        PB   = MatrixList(    zeros,   zeros,\
                              zeros,   ones)

        # U matrices, no units
        Uk   = MatrixList(  Ak/N.sqrt(2.)   ,  -Ak/N.sqrt(2.), \
                            ones/N.sqrt(2.) , ones/N.sqrt(2.) )

        Ukd  = MatrixList(    Ak_star/N.sqrt(2.),  ones/N.sqrt(2.),\
                             -Ak_star/N.sqrt(2.),  ones/N.sqrt(2.))

        Ukq  = MatrixList(  Akq/N.sqrt(2.)  ,  -Akq/N.sqrt(2.),\
                            ones/N.sqrt(2.) , ones/N.sqrt(2.) )

        Ukqd = MatrixList(   Akq_star/N.sqrt(2.),  ones/N.sqrt(2.),\
                            -Akq_star/N.sqrt(2.),  ones/N.sqrt(2.))

        # Q matrices - without units
        QA_mat = Ukqd*PA*Uk
        QB_mat = Ukqd*PB*Uk

        self.QA = QA_mat.return_list()
        self.QB = QB_mat.return_list()

        del(Id)
        del(SIG1)
        del(QA_mat)
        del(QB_mat)

        # Current matrix elements
        # units: e x hbar/(m a_0), fundamental unit of current
        current_prefactor = 2./3.*1j*hvF/Ha_to_eV 

        Hkx  = MatrixList(             zeros               ,  current_prefactor*Vkx,\
                      -current_prefactor*N.conjugate(Vkx)  ,          zeros         ) 

        Hky  = MatrixList(             zeros               ,  current_prefactor*Vky,\
                      -current_prefactor*N.conjugate(Vky)  ,          zeros         ) 

        Jx_mat = Ukd*Hkx*Uk
        Jy_mat = Ukd*Hky*Uk

        self.J_x = Jx_mat.return_list()
        self.J_y = Jy_mat.return_list()

        del(Jx_mat)
        del(Jy_mat)
        del(Hkx)
        del(Hky)

        # Electron-Phonon coupling
        # units: hbar^2/ m a_0^2, fundamental unit of energy

        elph_prefactor = N.sqrt(me_M)/N.sqrt(2.)*N.sqrt(Ha_to_eV/N.abs(self.hw_nu_q)) 

        # the u and v channels must be multiplied by (u x a0), in eV, to form the unitless
        # combination (u a0 / Ha).
        u = 1./Ha_to_eV 
        v = 1./Ha_to_eV 

        # E(-q) = E^*(q)
        E_star = N.conjugate(self.E_ph)

        m_11   = v*elph_prefactor* ( -E_star[3]*N.conjugate(Vqx)  \
                                     -E_star[4]*N.conjugate(Vqy)  )

        m_22   = v*elph_prefactor* (  E_star[0]*Vqx  \
                                     +E_star[1]*Vqy  )

        m_12   = u*elph_prefactor* (  E_star[0]*Vkqx   \
                                     +E_star[1]*Vkqy   \
                                     -E_star[3]*Vkx    \
                                     -E_star[4]*Vky  )

        m_21   = u*elph_prefactor* (  E_star[0]*N.conjugate(Vkx)   \
                                     +E_star[1]*N.conjugate(Vky)   \
                                     -E_star[3]*N.conjugate(Vkqx)  \
                                     -E_star[4]*N.conjugate(Vkqy)  )

        L_mat_u = MatrixList(  zeros,  m_12,\
                               m_21 ,  zeros)


        L_mat_v = MatrixList(  m_11 ,  zeros,
                               zeros,  m_22  )



        g_mat_u = Ukd*L_mat_u*Ukq
        g_mat_v = Ukd*L_mat_v*Ukq

        self.g_nu_u = g_mat_u.return_list()
        self.g_nu_v = g_mat_v.return_list()

        del(L_mat_u)
        del(L_mat_v)

        del(g_mat_u)
        del(g_mat_v)

        del(Uk)
        del(Ukd)
        del(Ukq)
        del(Ukqd)

        # Combining everything, in atomic units

        for alpha, J in zip(['x','y'],[self.J_x,self.J_y]):
            for L, QL in zip(['A','B'],[self.QA,self.QB]): 
                for n1 in [0,1]:
                    for n2 in [0,1]:
                        for n3 in [0,1]:
                            key   = (alpha,L, n1,n2,n3)
                            index = self.index_dictionary[key]

                            self.M[index,0,:] = J[n1,n2]*self.g_nu_u[n2,n3]*QL[n3,n1]
                            self.M[index,1,:] = J[n1,n2]*self.g_nu_v[n2,n3]*QL[n3,n1]



        return 
