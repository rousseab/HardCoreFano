#================================================================================
#
#                  Module Integrators
#           ===========================
#
#   The various classes implemented in this module are used to compute
#   integrations over the triangulation of the 1BZ.
#================================================================================

from module_Constants import *
from module_Grids import *

import sys


def AreaIntegrator(wedge,list_Fk):
    """
    This function integrates 
         I(w) = int d2k F_k(w) 

    assuming a linear dispersion within each triangle in the wedge.

    It is assumed that list_Fk = [ ---- Fk1 (w) ---- ]
                                 [ ---- Fk2 (w) ---- ]
                                 [ ---- Fk3 (w) ---- ]
                                 [ ---- ......  ---- ]
                                 [ ---- FkN (w) ---- ]
    """

    nk, nw = list_Fk.shape

    integral = complex(0.,0.)*N.zeros(nw)

    for tr, cross_product in zip(wedge.triangles_indices,wedge.cross_products):

        #F1 = list_Fk[tr[0]]
        F2 = list_Fk[tr[1]]
        F3 = list_Fk[tr[2]]

        integral += cross_product*(F2/3.+F3/6.)

    return integral
