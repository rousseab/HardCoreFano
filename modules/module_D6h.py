import numpy as N

def get_D6h():
    """
    This function codes explicitly the operations of the D6h group in 3 dimensions.
    This group leaves (for example) the benzene molecule invariant.
    The nomenclature below will follow the rules from Symmetry and Spectroscopy;
    in particular, see figure 1-25 of that book.

    The operations are:
         operations       class           description
    --------------------------------------------------------------------------------
        E               E              the identity

        C6, C6^{-1}    2C6             left and right rotation by pi/3 about the z axis

        C6^2, C6^{-2}  2C3             left and right rotation by 2 pi/3 about the z axis

        C6^3            C2             rotation by pi about the z axis

       C2'_{1,2,3}         3C2'            rotation by pi about the three axes which cross the center
                           of the hexagon and two corners.

       C2''_{1,2,3}        3C2''           rotation by pi about the three axes which cross the center
                           of the hexagon and two faces.

        I               I              the inversion

      sh*C6^2,sh*C6^{-2}   2S3             left and right rotation by 2 pi/3 about the z axis, followed by horizontal
                           reflection

      sh*C6,sh*C6^{-1}     2S6             left and right rotation by  pi/3 about the z axis, followed by horizontal
                           reflection

          sh               sh              reflection through horizontal plane

      sd_{1,2,3}           3sd             reflection through vertical planes crossing planes of hexagon

      sv_{1,2,3}           3sv             reflection through vertical planes crossing corners of hexagon


    """

    theta = N.pi/3
    c     = N.cos(theta)
    s     = N.sin(theta)

    # Identity
    E =  N.array([[  1.,  0.,  0.],              
                  [  0.,  1.,  0.],
                  [  0.,  0.,  1.]])


    # pi/3 rotation left
    C6_1  = N.array([[ c , s ,  0.], 
                     [-s , c ,  0.],   
                     [ 0., 0.,  1.]])

    # pi/3 rotation right
    C6_2  = N.array([[ c ,-s ,  0.], 
                     [ s , c ,  0.],   
                     [ 0., 0.,  1.]])


    # 2pi/3 rotation left
    C3_1  = N.dot(C6_1,C6_1)

    # 2pi/3 rotation right
    C3_2  = N.dot(C6_2,C6_2)

    # pi rotation 
    C2    = N.dot(C6_1,C3_1)

    # vertical mirror plane, and rotations by 2pi/3
    # these mirror planes cross the corners of the hexagon
    sv_1   =  N.array([[ -1.,  0.,  0.],              
                       [  0.,  1.,  0.], 
                       [  0.,  0.,  1.]])

    sv_2   =  N.dot(C3_1,sv_1)
    sv_3   =  N.dot(C3_2,sv_1)

    # mirror planes
    # these mirror planes cross the faces of the hexagon
    sd_1   =  N.dot(C6_1,sv_1)
    sd_2   =  N.dot(C3_1,sd_1)
    sd_3   =  N.dot(C3_2,sd_1)

    # reflection about the horizontal plane
    sh     =  N.array([[  1.,  0.,  0.],              
                       [  0.,  1.,  0.], 
                       [  0.,  0., -1.]])

    # Inversion
    I      =  N.array([[ -1.,  0.,  0.],              
                       [  0., -1.,  0.], 
                       [  0.,  0., -1.]])


    # rotation by pi about the y axis, which crosses the corners,
    # followed by 2pi/3 rotation of the C2 axis
    C2p_1  =  N.array([[ -1.,  0.,  0.],              
                       [  0.,  1.,  0.], 
                       [  0.,  0., -1.]])


    C2p_2  =  N.dot(C3_1,C2p_1)
    C2p_3  =  N.dot(C3_2,C2p_1)


    # rotation by pi about the axes crossing the faces of the hexagon
    C2pp_1 = N.dot(C6_1,C2p_1)
    C2pp_2 = N.dot(C3_1,C2pp_1)
    C2pp_3 = N.dot(C3_2,C2pp_1)


    # improper rotations
    S3_1   = N.dot(sh,C3_1)
    S3_2   = N.dot(sh,C3_2)

    S6_1   = N.dot(sh,C6_1)
    S6_2   = N.dot(sh,C6_2)


    list_D6h = N.array([ E, C6_1, C6_2, C3_1, C3_2, C2, C2p_1, C2p_2, C2p_3, C2pp_1, C2pp_2, C2pp_3,
                 I, S3_1, S3_2, S6_1, S6_2, sh,  sd_1,  sd_2,  sd_3,   sv_1,   sv_2,   sv_3]) 


    return list_D6h

def get_D6h_character_of_representations():

    A1g = N.array([ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ])
    A2g = N.array([ 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1 ])
    B1g = N.array([ 1,-1,-1, 1, 1,-1, 1, 1, 1,-1,-1,-1, 1,-1,-1, 1, 1,-1, 1, 1, 1,-1,-1,-1 ])
    B2g = N.array([ 1,-1,-1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1, 1, 1 ])
    E1g = N.array([ 2, 1, 1,-1,-1,-2, 0, 0, 0, 0, 0, 0, 2, 1, 1,-1,-1,-2, 0, 0, 0, 0, 0, 0 ])
    E2g = N.array([ 2,-1,-1,-1,-1, 2, 0, 0, 0, 0, 0, 0, 2,-1,-1,-1,-1, 2, 0, 0, 0, 0, 0, 0 ])

    A1u = N.array([ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 ])
    A2u = N.array([ 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1 ])
    B1u = N.array([ 1,-1,-1, 1, 1,-1, 1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1, 1,-1,-1,-1, 1, 1, 1 ])
    B2u = N.array([ 1,-1,-1, 1, 1,-1,-1,-1,-1, 1, 1, 1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1 ])
    E1u = N.array([ 2, 1, 1,-1,-1,-2, 0, 0, 0, 0, 0, 0,-2,-1,-1, 1, 1, 2, 0, 0, 0, 0, 0, 0 ])
    E2u = N.array([ 2,-1,-1,-1,-1, 2, 0, 0, 0, 0, 0, 0,-2, 1, 1, 1, 1,-2, 0, 0, 0, 0, 0, 0 ])



    representations = {'A1g':A1g,'A2g':A2g,'B1g':B1g,'B2g':B2g,'E1g':E1g,'E2g':E2g,
               'A1u':A1u,'A2u':A2u,'B1u':B1u,'B2u':B2u,'E1u':E1u,'E2u':E2u}

    return representations 


def check_D6h():

    list_D6h = get_D6h()

    NG = 1.*len(list_D6h)

    traces  = list_D6h[:,0,0]+list_D6h[:,1,1]+list_D6h[:,2,2]

    representations = get_D6h_character_of_representations()


    NR = 1.*len(representations)

    print '\nTest 1: 3D Representation'
    print '\tGroup matrices:'
    for rep in representations:

        G = representations[rep]
        o = N.dot(G,traces)/NG

        print '\trepresentation %s, dot product = %+4.3f'%(rep,o)

    print '\nTest 2: Orthonormality of group representations'


    M = N.zeros([NR,NR])
    for i,r1 in enumerate(representations):
        G1 = representations[r1]
        for j,r2 in enumerate(representations):
            G2 = representations[r2]
            M[i,j] = N.dot(G1,G2)/NG

    print '\t product matrix :'
    print M


    print '\t| product matrix - identity | = %12.8f'%(N.linalg.norm(M-N.eye(NR)))

    return



#check_D6h()

