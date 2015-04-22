import common
reload(common)
from common import *


print "================================================================================"
print "#                                                                               "
print "#    This test computes                                                         "
print "#                                                                               "
print "#  I^{alpha,nu}_{n1n2n3n4} = sum_{k} J^{alpha}_{n1n2}(k) g^{nu}_{n3n4}(k,k)     "
print "#                                                                               "
print "#  The Kuzmenko susceptibility takes the form                                   "
print "#                                                                               "
print "#   chi^{alpha,nu} ~ sum_{k} sum_{n1 n2}                                        "
print "#             J^{alpha}_{n1n2}(k) g^{nu}_{n2 n1}(k) G_{n1}(k) G_{n2}(k)         "
print "#                                                                               "
print "#  where the Green functions G(k) are diagonal and can be assumed to have the   "
print "#  full symmetry of the k-lattice. diagram vanishes in graphene by symmetry.    "
print "#                                                                               "
print "# Thus, if I^{alpha,nu}_{n1n2n2n1} = 0, chi^{alpha,nu} must vanish.             "
print "================================================================================"

path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_10.0/modules/'
list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)

q_vector = N.array([0.,0.])

list_hw, list_hw_renormalized, Eq = Renormalizor.get_frequency_and_polarization(q_vector)


print 'q  = ( %8.4f,  %8.4f)'%(q_vector[0],q_vector[1])


print "# ====  Computing  ===="
print "#  I^{alpha,nu}_{n1n2n3n4}" 

result_dic = {}
for alpha in ["x","y"]:
    for nu in N.arange(6):
        for n1 in [0,1]:
            for n2 in [0,1]:
                for n3 in [0,1]:
                    for n4 in [0,1]:

                        key = (alpha,nu,n1,n2,n3,n4)
                        result_dic[key] = 0.


nmax_coarse = 16
nmax_fine   = 64
n_blocks_coarse_to_fine = 4
include_Gamma = True

grid = TesselationDoubleGrid(nmax_coarse, nmax_fine, n_blocks_coarse_to_fine,include_Gamma, clip_grid=False)


for nu in N.arange(0,6):
    print 'nu = %i'%nu
    E_phonon_polarization = Eq[:,nu]
    hw_nu_q               = list_hw[nu]

    for i,wedge in enumerate(grid.list_wedges):

        M = MGridFunction(q_vector,E_phonon_polarization ,hw_nu_q, wedge)

        for alpha, J_alpha in zip(["x","y"],[M.J_x,M.J_y]):
            for n1 in [0,1]:
                for n2 in [0,1]:
                    for n3 in [0,1]:
                        for n4 in [0,1]:

                            key = (alpha,nu,n1,n2,n3,n4)
                            result_dic[key] += N.sum(J_alpha[n1,n2,:]*M.g_nu_u[n3,n4,:])


# compute the traces

tol = 1e-12

print "# nu  alpha   n1   n2   |I^{alpha,nu}_{n1n2n2n1}|"
print "#================================================"

for nu in N.arange(0,6):
    for alpha in ["x","y"]:

            for n1 in [0,1]:
                for n2 in [0,1]:
                    key = (alpha,nu,n1,n2,n2,n1)
                    T = N.abs(result_dic[key])

                    print "  %i     %s     %i    %i      %10.3e "%(nu,alpha,n1,n2,T)
                    if T > tol:
                        print '**** ERROR! THE TRACE IS NON VANISHING! ****'
                        sys.exit()


print "\n\n  Trace is always numerically zero: test is passed"
