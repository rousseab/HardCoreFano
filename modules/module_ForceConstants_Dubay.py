#================================================================================
#
#           Module Force Constants Dubay
#           ===========================
#
#   This module implements the force constants to extract the phonon
#   dispersions and polarizations, following the results of Dubay and Kresse.
#
#
#================================================================================
from module_Constants import *
from module_D6h import *

import socket

def getDubayForceConstants():

    hostname = socket.gethostname()

    if hostname == 'ferron.pmc.umontreal.ca':
        datafile = open('/Users/Shared/Bruno/work/projects/Graphene_Fano/python_work/HardCoreKernelFano/modules/force_constants.txt')

    else:
        datafile = open('/RQusagers/roussea4/python_depository/HCF_5.0/modules/force_constants.txt')


    DubayForceConstants_cartesian = {}
    DubayForceConstants_radial    = {}

    for line in datafile:

        if line[0] != '#':

            strip_line = line.strip().split()

            label = '%3s'%(strip_line[0].replace('\'','p'))

            dat = []
            for el in strip_line[1:]:
                dat.append(-float(el)*bohr_in_angst**2)
        
            Matrix_cartesian = N.array([[dat[0],dat[2],   0.],
                                        [dat[1],dat[3],   0.],
                                        [    0.,    0.,dat[8]]])

            Matrix_radial    = N.array([[dat[4],dat[6],   0.],
                                        [dat[5],dat[7],   0.],
                                        [    0.,    0.,dat[8]]])

            DubayForceConstants_cartesian[label] = Matrix_cartesian 
            DubayForceConstants_radial   [label] = Matrix_radial


    return DubayForceConstants_cartesian, DubayForceConstants_radial   

def build_real_space_lattice(nmax):
    """
    Builds an array of real space vectors, R, in a box [-nmax, nax],[-nmax,nmax].
    This is overkill as it generates more than what we need, but it is simple
    and robust.
    """

    list_i = N.arange(-nmax,nmax+1)

    I, J   = N.meshgrid(list_i,list_i)

    If = I.flatten()
    Jf = J.flatten()

    list_R = If[:,N.newaxis]*a1+Jf[:,N.newaxis]*a2


    return list_R

def create_neigbhors_in_wedge(list_carbon_positions):
    """
    Extract only the positions inside a given wedge, between 30 and 90 degrees, 
    in order to reproduce the Dubay tree of points. 
    """

    tol = 1e-8
    list_p = list_carbon_positions-tau_B

    list_p_wedge = []

    for p in list_p:

        append = False

        if p[0] >= 0.:

            if N.abs(p[0]) < tol:

                if p[1] >= 0.:
                    append = True

            else:
                theta = N.arctan(p[1]/p[0])

                if theta-N.pi/6. >= -tol  and theta <= N.pi/2.:
                    append = True

        if append:
            list_p_wedge.append(p)


    list_p_wedge = N.array(list_p_wedge)

    norm = N.sum(list_p_wedge**2 ,axis=1)
    I    = N.argsort(norm)

    return list_p_wedge[I]

def get_D3():
    """
    This function codes explicitly the operations of the D3 group in 2 dimensions.


    The operations are:
         operations       class           description
    --------------------------------------------------------------------------------
        E               E              the identity

      C3_1, C3_2           2C3             left and right rotation by 2 pi/3 about the z axis

      sv_{1,2,3}           3sv             reflection through vertical planes 


    """

    theta = 2.*N.pi/3
    c     = N.cos(theta)
    s     = N.sin(theta)

    # Identity
    E =  N.array([[  1.,  0.],              
              [  0.,  1.]])

    # 2pi/3 rotation left
    C3_1  = N.array([[ c , s ], 
             [-s , c ]])

    # pi/3 rotation right
    C3_2  = N.array([[ c ,-s ], 
             [ s , c ]])


    # vertical mirror plane, and rotations by 2pi/3
    sv_1   =  N.array([[ -1.,  0.],
               [  0.,  1.]])

    sv_2   =  N.dot(C3_1,sv_1)
    sv_3   =  N.dot(C3_2,sv_1)


    list_D3 = N.array([ E, C3_1, C3_2, sv_1, sv_2, sv_3]) 


    return list_D3

def build_force_constants_Dubay():

    nmax = 5
    list_R = build_real_space_lattice(nmax)

    DFC_cartesian, DFC_radial = getDubayForceConstants()

    tol = 1e-8


    list_RA = list_R+tau_A
    list_RB = list_R+tau_B

    list_C  = N.concatenate([list_RA,list_RB])

    nmax = 26
    list_p_wedge = create_neigbhors_in_wedge(list_C)[1:nmax]
    list_labels  = ['  1','  2','  3','  4','  5',' 5p','  6' ,'  7','  8','  9','10p',' 10',' 11',
                    ' 12',' 13',' 14','15p',' 15',' 16','17p',' 17',' 18',' 19','20p',' 20']


    list_D3 = get_D3()

    list_total_labels = []
    list_p            = []
    list_sigma        = []

    list_S2 = []
    for label, p_wedge in zip(list_labels,list_p_wedge):

        # Action of symmetry elements on position p_wedge
        Dp = N.dot(list_D3,p_wedge)

        # Compute all norms || Dp_i - Dp_j ||^2
        x  = Dp[:,0]
        y  = Dp[:,1]

        Ax,Bx = N.meshgrid(x,x)
        Ay,By = N.meshgrid(y,y)

        # Add lower triangular matrix, to avoid zero norm below the diagonal
        S2 = (Ax-Bx)**2+(Ay-By)**2+N.tri(6,6,0)


        # compute array describing which pairs of vectors are equal, if any
        Indices_equal = N.argwhere(S2 < tol)

        # find the indices of redundant vectors, 
        indices_redundant = Indices_equal[:,1]

        truth = N.array([True,True,True,True,True,True])

        truth[indices_redundant] = False

        for dp,sigma in zip(Dp[truth,:],list_D3[truth]):
            list_p.append(dp)
            list_sigma.append(sigma)
            list_total_labels.append(label)


    SIGMA = N.array(list_sigma)

    SB  = N.array([[1.,0., 0.],[0.,1.,0.],[0.,0.,1.]])
    PB  = N.array(list_p)
    L   = list_total_labels

    PA = deepcopy(PB)

    SA  = N.array([[1.,0., 0.],[0.,-1.,0.],[0.,0.,1.]])
    PA[:,1] = -PA[:,1]

    # Create the force constants!
    dic_DR = {}

    for P, S, kappa1, tau_kappa1 in zip([PA,PB],[SA,SB],[0,1],[tau_A,tau_B]):

        for p, sigma, lab in zip(P,SIGMA, L):

            # use cartesian force constants, as the radial ones
            # appears to be ambiguously defined
            FC0 = DFC_cartesian[lab]

            # space transformation to go from Dubay's configuration to this point
            sig3D = N.diag([1.,1.,1.])
            sig3D[:2,:2] = sigma

            U     = N.dot(S,sig3D)

            #FC  = N.dot(N.transpose(U),N.dot(FC0,U))
            FC  = N.dot(U,N.dot(FC0,N.transpose(U)))

            # for each distance p, find out R2, kappa2 

            # p = [R2 + tau_kappa2 ] -tau_kappa1

            d = p+tau_kappa1

            RA = d-tau_A
            RB = d-tau_B

            i1A = N.dot(b1,RA)/(2.*N.pi)
            i2A = N.dot(b2,RA)/(2.*N.pi)

            i1B = N.dot(b1,RB)/(2.*N.pi)
            i2B = N.dot(b2,RB)/(2.*N.pi)

            if N.sqrt((i1A - N.round(i1A))**2 +  (i2A - N.round(i2A))**2) < tol:
                i1 = N.int(N.round(i1A))
                i2 = N.int(N.round(i2A))
                kappa2 = 0
            elif N.sqrt((i1B - N.round(i1B))**2 +  (i2B - N.round(i2B))**2) < tol:
                i1 = N.int(N.round(i1B))
                i2 = N.int(N.round(i2B))
                kappa2 = 1

            else:
                print 'PROBLEM!'


            DR = N.zeros([6,6])
            if kappa1 == 0 and kappa2 == 0:
                DR[:3,:3] += FC
            elif kappa1 == 0 and kappa2 == 1:
                DR[:3,3:] += FC
            elif kappa1 == 1 and kappa2 == 0:
                DR[3:,:3] += FC
            elif kappa1 == 1 and kappa2 == 1:
                DR[3:,3:] += FC


            # contribution D
            #indices_R1_minus_R2 = (-i1,-i2)
            indices_R1_minus_R2 = (-i1,-i2)


            if dic_DR.has_key(indices_R1_minus_R2):
                dic_DR[indices_R1_minus_R2] += DR
            else:
                dic_DR[indices_R1_minus_R2]  = DR


    # Extract explicit arrays
    m = carbon_mass/electron_mass

    list_R  = []
    list_FC = []
    for indices_R1_minus_R2, DR in dic_DR.iteritems():

        R = indices_R1_minus_R2[0]*a1+indices_R1_minus_R2[1]*a2
        list_R.append(R)

        FC = DR/m
        list_FC.append(FC)


    list_R = N.array(list_R)
    I      = N.argsort(N.sum(list_R**2,axis=1))
    list_R = list_R[I]
    list_FC= N.array(list_FC)[I]

    # impose sum rule
    SR1 = -N.sum(list_FC[:,:3,:3]+ list_FC[:,:3,3:],axis=0)
    SR2 = -N.sum(list_FC[:,3:,:3]+ list_FC[:,3:,3:],axis=0)

    list_FC[0,:3,:3] = SR1
    list_FC[0,3:,3:] = SR2


    return list_R, list_FC

def build_force_constants_Dubay_using_symmetry():


    list_D6h = get_D6h()

    NG = len(list_D6h)


    tol = 1e-12

    #================================================================================
    # get the minimal set of force constants
    #================================================================================

    nmax = 5
    list_R = build_real_space_lattice(nmax)

    DFC_cartesian, DFC_radial = getDubayForceConstants()


    list_RA = list_R+tau_A
    list_RB = list_R+tau_B

    list_C  = N.concatenate([list_RA,list_RB])

    nmax = 26
    list_d_wedge = create_neigbhors_in_wedge(list_C)[1:nmax]
    # labels, ordered in the right way given the sorting done in the routine before
    list_labels  = ['  1','  2','  3','  4','  5',' 5p','  6' ,'  7','  8','  9','10p',' 10',' 11',
            ' 12',' 13',' 14','15p',' 15',' 16','17p',' 17',' 18',' 19','20p',' 20']


    # Put the force constants in a dictionary
    dic_D0 = {}

    kappa1 = 1 # the force constants are given for tau_B at the origin

    for d, label in zip(list_d_wedge,list_labels):

        # define a key which contains R1-R2, kappa1, kappa2

        r = d+tau_B

        rA = r-tau_A
        rB = r-tau_B

        ib1A = N.dot(rA,b1)/(2.*N.pi)
        ib2A = N.dot(rA,b2)/(2.*N.pi)

        testA = N.sqrt( (ib1A-N.round(ib1A))**2 + (ib2A-N.round(ib2A))**2 )  < tol

        ib1B = N.dot(rB,b1)/(2.*N.pi)
        ib2B = N.dot(rB,b2)/(2.*N.pi)

        testB = N.sqrt( (ib1B-N.round(ib1B))**2 + (ib2B-N.round(ib2B))**2 )  < tol

        if testA:
            ib1    = N.int(N.round(ib1A))
            ib2    = N.int(N.round(ib2A))
            kappa2 = 0
        elif testB:
            ib1    = N.int(N.round(ib1B))
            ib2    = N.int(N.round(ib2B))
            kappa2 = 1

        else:
            print 'ERROR!'


        key = ( -ib1, -ib2, kappa1, kappa2 ) # R1-R2, kappa1, kappa2


        if dic_D0.has_key(key):
            print 'Dictionary dic_D0 already has key!'

        dic_D0[key] = [ DFC_cartesian[label], label ] 


    #================================================================================
    # Generate all the permutation matrices
    #================================================================================

    list_Pg = []

    for g in list_D6h:

        Rg = g[:2,:2]
        Pg = N.zeros([2,2])
        for kappa1, tau_kappa1 in enumerate([tau_A,tau_B]):
            for kappa2, tau_kappa2 in enumerate([tau_A,tau_B]):

                R = tau_kappa1 - N.dot(Rg,tau_kappa2)

                ib1 = N.dot(R,b1)/(2.*N.pi)
                ib2 = N.dot(R,b2)/(2.*N.pi)

                error = (ib1 - N.round(ib1))**2 + (ib2 - N.round(ib2))**2

                if error < tol:
                    Pg[kappa1,kappa2] = 1 

        list_Pg.append(Pg)

    list_Pg = N.array(list_Pg )

        #================================================================================
        # Fill force constant dictionary with all possibilities
        #================================================================================


    # Sort all contributions in sets of repetitions
    dic_repetitions = {}
    for key0, item0 in dic_D0.iteritems():

        D0    = item0[0]
        label = item0[1]

        R      = key0[0]*a1 + key0[1]*a2
        kappa1 = key0[2]
        kappa2 = key0[3]

        for Rg_3D,Pg in zip(list_D6h,list_Pg):

            Rg = Rg_3D[:2,:2]

            D  = N.dot(Rg_3D, N.dot(D0,Rg_3D.transpose()))

            # figure out the basis indices
            if Pg[0,kappa1] == 1:
                bar_kappa1 = 0
                tau_bar_kappa1 = tau_A

            elif Pg[1,kappa1] == 1:
                bar_kappa1 = 1
                tau_bar_kappa1 = tau_B
            else:
                print 'PROBLEM!!!'

            if Pg[0,kappa2] == 1:
                bar_kappa2 = 0
                tau_bar_kappa2 = tau_A

            elif Pg[1,kappa2] == 1:
                bar_kappa2 = 1
                tau_bar_kappa2 = tau_B
            else:
                print 'PROBLEM!!!'


            # figure out the new R vector


            RgA = N.dot(Rg,tau_A)
            RgB = N.dot(Rg,tau_B)
                
            delta_R = (Pg[bar_kappa1,0]-Pg[bar_kappa2,0])*RgA+(Pg[bar_kappa1,1]-Pg[bar_kappa2,1])*RgB \
                    - (tau_bar_kappa1-tau_bar_kappa2)


            RgR = N.dot(Rg,R)
            Rp  = RgR +delta_R 


            # Build new key
            ib1_f = N.dot(Rp,b1)/(2.*N.pi)
            ib2_f = N.dot(Rp,b2)/(2.*N.pi)

            ib1 = N.int(N.round(N.dot(Rp,b1)/(2.*N.pi)))
            ib2 = N.int(N.round(N.dot(Rp,b2)/(2.*N.pi)))

            error = N.sqrt((ib1-ib1_f)**2+(ib2-ib2_f)**2)

            if error > tol:
                print 'PROBLEM!'
                bloup

            key = ( ib1, ib2, bar_kappa1, bar_kappa2 ) # R1-R2, kappa1, kappa2

            bundle = [R,kappa1,kappa2,Rg_3D, Pg,RgR, delta_R, Rp, bar_kappa1, bar_kappa2, D]

            if dic_repetitions.has_key(key):
                dic_repetitions[key].append(bundle)
            else:
                dic_repetitions[key] = [bundle]


    debug = False
    if debug:
        for key, item in dic_repetitions.iteritems():
            print '%i contributions to same point! '%len(item)
            for bundle in item:
                # bundle = [R,kappa1,kappa2,Rg_3D, Pg,RgR, delta_R, Rp, bar_kappa1,bar_kappa2]
                R,kappa1,kappa2,Rg_3D, Pg,RgR, delta_R, Rp, bar_kappa1,bar_kappa2, D = bundle 
                print '\t\t symmetry operation: '
                print '\t\t  R(g) = [ %+8.4f  %+8.4f   %+8.4f ] '%(Rg_3D[0,0],Rg_3D[0,1],Rg_3D[0,2])
                print '\t\t         [ %+8.4f  %+8.4f   %+8.4f ] '%(Rg_3D[1,0],Rg_3D[1,1],Rg_3D[1,2])
                print '\t\t         [ %+8.4f  %+8.4f   %+8.4f ] '%(Rg_3D[2,0],Rg_3D[2,1],Rg_3D[2,2])
                print '\n'
                print '\t\t     R  = [ %+8.4f %+8.4f ] '%(R[0],R[1])
                print '\t\t R(g).R = [ %+8.4f %+8.4f ] '%(RgR[0],RgR[1])
                print '\t\t  Del_R = [ %+8.4f %+8.4f ] '%(delta_R[0],delta_R[1])
                print '\t\t     Rp = [ %+8.4f %+8.4f ] '%(Rp[0],Rp[1])
                print '\t\t kappa1 = %i   kappa2 = %i --> bar_kappa1 = %i  bar_kappa2 = %i'%(kappa1,kappa2,bar_kappa1,bar_kappa2)
                print '\n'

    # Put average results in a dictionary; the force constants are only stated to 4 decimals; there may be errors
    # that need to be whipped out by averaging
    dic_D = {}

    for key, item in dic_repetitions.iteritems():
        D = N.zeros([3,3])
        for bundle in item:
            D += bundle[-1]/len(item)

        dic_D[key] = D


    # Build the array of force constants
    list_R = []
    for key in dic_D.keys():
        R = key[0]*a1 + key[1]*a2

        Found = False
        for Rp in list_R:
            if N.linalg.norm(R-Rp) < tol:
                Found = True

        if not Found:
            list_R.append(R)

    list_R = N.array(list_R)
    NR     = len(list_R)

    norm2  = N.sum(list_R**2,axis=1)
    I      = N.argsort(norm2)
    list_R = list_R[I]



    list_FC = N.zeros([NR,6,6])

    m = carbon_mass/electron_mass

    for key, D in dic_D.iteritems():

        R = key[0]*a1 + key[1]*a2

        kappa1 = key[2]
        kappa2 = key[3]

        key_transpose = (-key[0],-key[1],kappa2,kappa1)

        D_transpose   = N.transpose(dic_D[key_transpose])

        # find the index of R

        for i, Rp in enumerate(list_R):

            if N.linalg.norm(R-Rp) < tol:
                break

        # Impose transpose symmetry
        FC = 0.5*(D+D_transpose)/m

        if kappa1 == 0 and kappa2 == 0:
            list_FC[i,:3,:3] = FC

        if kappa1 == 0 and kappa2 == 1:
            list_FC[i,:3,3:] = FC

        if kappa1 == 1 and kappa2 == 0:
            list_FC[i,3:,:3] = FC
     
        if kappa1 == 1 and kappa2 == 1:
            list_FC[i,3:,3:] = FC
         

    # impose sum rule
    SR1 = -N.sum(list_FC[:,:3,:3]+ list_FC[:,:3,3:],axis=0)
    SR2 = -N.sum(list_FC[:,3:,:3]+ list_FC[:,3:,3:],axis=0)

    list_FC[0,:3,:3] = SR1
    list_FC[0,3:,3:] = SR2


    return list_R, list_FC
