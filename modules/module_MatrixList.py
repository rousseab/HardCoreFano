#================================================================================
#
#           Module MatrixList
#           =================
#
#       This module implements an object which can store a list of matrices
#       and operate on them in an intuitive manner.
#
#================================================================================

from module_Constants import *

class MatrixList:
    """
    The purpose of this class is to store and manipulate
    lists of 2x2 matrices. This will allow the direct
    implementation of equations in the graphene project.
    """

    def __init__(self,m_11,m_12,m_21,m_22):
        """
        Initialize the object with components
            [ m_11  m_12 ]
            [ m_21  m_22 ]
        where each component is assumed to be a list of N elements.
        """
        self.m_11 = deepcopy(m_11)
        self.m_12 = deepcopy(m_12)
        self.m_21 = deepcopy(m_21)
        self.m_22 = deepcopy(m_22)

    def return_list(self):

        return N.array([ [self.m_11, self.m_12],[self.m_21, self.m_22]])

    def __mul__(self,B_matrix):

        new_11 = self.m_11*B_matrix.m_11 + self.m_12*B_matrix.m_21 

        new_12 = self.m_11*B_matrix.m_12 + self.m_12*B_matrix.m_22 

        new_21 = self.m_21*B_matrix.m_11 + self.m_22*B_matrix.m_21 

        new_22 = self.m_21*B_matrix.m_12 + self.m_22*B_matrix.m_22

        New = MatrixList(new_11,new_12,new_21, new_22)

        return New

