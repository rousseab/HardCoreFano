from common import *



def get_frequencies_and_polarizations_LOCAL(q,list_FC,list_R):

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

    return list_hw, Eq, Bq

def fix_phase(vector):

    tol = 1e-12
    for x in vector:
        if N.abs(x) > tol:
            phase = x/N.abs(x)
            break

    new_vector = vector*N.conjugate(phase)
    return new_vector 



list_R, list_FC  =  build_force_constants_Dubay_using_symmetry()

nkpath = 500

kpath = Path(acell,nkpath)

Kplus  = kpath.K
Kminus = -kpath.K
G      = kpath.G
M      = kpath.M

list_S = [G,M,Kplus,Kminus]

path='/Users/Bruno/work/Projects/fano_project/HardCoreKernelFano_3.0/modules/'
Renormalizor = CompteMauriRenormalize(path,list_FC,list_R)

q = Kplus

phaseA = N.exp(-1j*N.dot(q,tau_A))
phaseB = N.exp(-1j*N.dot(q,tau_B))

list_hw_LOC, Eq_LOC, Bq_LOC = get_frequencies_and_polarizations_LOCAL(q,list_FC,list_R)
list_hw, list_hw_renormalized, Eq =  Renormalizor.get_frequency_and_polarization(q)

bq = []

for v in N.transpose(Bq_LOC):
    nv = fix_phase(v)

    print N.round(nv,4)


