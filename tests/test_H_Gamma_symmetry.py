#================================================================================
# Plot a nice tesselation of the 1BZ
#
#================================================================================
import common
reload(common)
from common import *


print "================================================================================"
print "#                                                                               "
print "#    This test computes                                                         "
print "#                                                                               "
print "#  I^{alpha,nu,L}_{n1n2n3n4n5n6} =                                              "
print "#        sum_{k} J^{alpha}_{n1n2}(k) g^{nu}_{n3n4}(k,k) Q^L_{n5n6}(k,k)         "
print "#                                                                               "
print "#  The H function at Gamma will depend on I, such that                          "
print "#   H^{alpha,nu,L}(Gamma) ~ sum_{k} sum_{n1 n2 n3}                              "
print "#          J^{alpha}_{n1n2}(k) g^{nu}_{n2 n3}(k)  Q^L_{n3n1} x (Green functions)"
print "#                                                                               "
print "#  where the Green functions G(k) are diagonal and can be assumed to have the   "
print "#  full symmetry of the k-lattice.                                              "
print "#                                                                               "
print "# Thus, if I^{alpha,nu}_{n1n2n2n2n3n3n1} != 0, H^{alpha,nuL}(Gamma) != 0.       "
print "================================================================================"


path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_10.0/modules/'
list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)

q_vector = N.array([0.,0.])

list_hw, list_hw_renormalized, Eq = Renormalizor.get_frequency_and_polarization(q_vector)


print 'q  = ( %8.4f,  %8.4f)'%(q_vector[0],q_vector[1])


print "# ====  Computing  ===="
print "#"
print "#  I^{alpha,nu,L}_{n1n2n3n4n5n6} = sum_{k} J^{alpha}_{n1n2}(k) g^{nu}_{n3n4}(k,k) Q^L_{n5n6}(k,k)"
print "#"
print "#====================="

result_dic = {}
for alpha in ["x","y"]:
    for L in ["A","B"]:
        for nu in N.arange(6):
            for n1 in [0,1]:
                for n2 in [0,1]:
                    for n3 in [0,1]:
                        for n4 in [0,1]:
                            for n5 in [0,1]:
                                for n6 in [0,1]:

                                    key = (alpha,L,nu,n1,n2,n3,n4,n5,n6)
                                    result_dic[key] = 0.

nmax =8 
grid = TesselationGrid(nmax)

for nu in N.arange(0,6):
    print 'nu = %i'%nu
    E_phonon_polarization = Eq[:,nu]
    hw_nu_q               = list_hw[nu]


    for i,wedge in enumerate(grid.list_wedges):

        M = MGridFunction(q_vector,E_phonon_polarization ,hw_nu_q, wedge)

        for alpha, J_alpha in zip(["x","y"],[M.J_x,M.J_y]):
            for L, QL in zip(["A","B"],[M.QA,M.QB]):

                for n1 in [0,1]:
                    for n2 in [0,1]:
                        for n3 in [0,1]:
                            for n4 in [0,1]:
                                for n5 in [0,1]:
                                    for n6 in [0,1]:

                                        key = (alpha,L,nu,n1,n2,n3,n4,n5,n6)

                                        result_dic[key] += N.sum(J_alpha[n1,n2,:]*M.g_nu_u[n3,n4,:]*QL[n5,n6,:])


# compute the trace elements

tol = 1e-12

print  "NON VANISHING ELEMENTS:" 
print "# nu  alpha   L   n1  n2  n3   ABS[I] "
print "#======================================"



for nu in N.arange(0,6):
    for alpha in ["x","y"]:
        for L in ["A","B"]:

            for n1 in [0,1]:
                for n2 in [0,1]:
                    for n3 in [0,1]:

                        key = (alpha,L,nu,n1,n2,n2,n3,n3,n1)
                        T   = N.abs(result_dic[key])

                        if T > tol:
                            print "  %i    %s       %s   %i   %i  %i  %10.3e"%(nu,alpha,L,n1,n2,n3,T)

