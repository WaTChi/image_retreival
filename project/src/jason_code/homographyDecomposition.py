# a lot of these formulas are from:
# Motion and Structure From Motion in a Piecewise Planar Environment
# by O. D. Faugeras, F. Lustman

import itertools
import math
import cv
import numpy as np

def extract_singular_values(matrix):
    result = []
    for i in range(min(matrix.shape)):
        result.append(matrix[i,i])
    return result

def num_unique(ls):
    len(set(ls))

def three_distinct(D):
    # first, the positive case
    # there are two values that can both be +/- 1, so use all combinations
    e_vals = [i for i in itertools.product((-1, 1), (-1, 1))]
    D_sq = [val * val for val in D]
    norm_candidates = \
      [(E[0]*math.sqrt((D_sq[0]-D_sq[1])/(D_sq[0]-D_sq[2])),
        0,
        E[1]*math.sqrt((D_sq[1]-D_sq[2])/(D_sq[0]-D_sq[2]))) for E in e_vals]
    trig_candidates = \
      [(E[0]*E[1]*math.sqrt((D_sq[0]-D_sq[1])*(D_sq[1]-D_sq[2]))/((D[0]+D[2])*D[1]),
        (D_sq[1]+D[0]*D[2])/((D[0]+D[2])*D[1])) for E in e_vals]
    R_candidates = \
      [np.mat([[trig[1], 0, -trig[0]],
               [0, 1, 0],
               [trig[0], 0, trig[1]]]) for trig in trig_candidates]
    T_candidates = \
      [(D[0]-D[2])*np.mat([[norm[0]],[0],[-norm[2]]]) for norm in norm_candidates]
    # and now, the negative case
    trig_candidates2 = \
      [(E[0]*E[1]*math.sqrt((D_sq[0]-D_sq[1])*(D_sq[1]-D_sq[2]))/((D[0]-D[2])*D[1]),
        (D[0]*D[2]-D_sq[1])/((D[0]-D[2])*D[1])) for E in e_vals]
    R_candidates2 = \
      [np.mat([[trig[1], 0, trig[0]],
               [0, 1, 0],
               [trig[0], 0, -trig[1]]]) for trig in trig_candidates2]
    T_candidates2 = \
      [(D[0]+D[2])*np.mat([[norm[0]],[0],[norm[2]]]) for norm in norm_candidates]
    return R_candidates+R_candidates2, T_candidates+T_candidates2

def find_SVD(homography):
    W = cv.CreateMat(3, 3, cv.CV_64F)
    U = cv.CreateMat(3, 3, cv.CV_64F)
    V = cv.CreateMat(3, 3, cv.CV_64F)
    cv.SVD(homography, W, U, V, cv.CV_SVD_U_T | cv.CV_SVD_V_T)
    return W, U, V
    
def euler_angles(R):
    return 180*math.acos(R[2,2])/math.pi, 180*math.atan2(R[2,0],R[2,1])/math.pi, 180*-math.atan2(R[0,2],R[1,2])/math.pi
    #return 180*math.acos(R[2,2])/math.pi, 180*math.atan2(R[0,2],R[1,2])/math.pi, 180*-math.atan2(R[2,0],R[2,1])/math.pi

def extract_R_and_T(homography):
    A, U, V = [np.asarray(mat) for mat in find_SVD(homography)]
    s = np.linalg.det(U) * np.linalg.det(V)
    D = extract_singular_values(A)
    results = []
    rotations, translations = three_distinct(D)
    #print "----"
    #print rotations
    #print "----"
    for values in zip(rotations, translations):
        #print values[0]
        #print "------"
        #results.append((s*U*values[0]*np.transpose(V), U*values[1]))
        results.append((s*U*values[0]*V, U*values[1]))
    for result in results:
        print "-----"
        print "Rotation: (pitch, yaw, roll) "
        print euler_angles(result[0])
        print "Translation: "
        print result[1]

  
