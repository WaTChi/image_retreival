#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="aaronh"
__date__ ="$Aug 29, 2011 5:02:12 PM$"

import geom
import numpy as np
import numpy.linalg as alg
import numpy.random as rand
import util
import cv
import homographyDecomposition
import computePose

if __name__ == "__main__":

    # camera 1
    w1, h1, fov1 = 768, 512, 60*np.pi/180
    K1 = geom.cameramat(w1,h1,fov1)
    #y1, p1, r1 = 75, -2, -1
    #R1 = geom.Rfrom_ypr(y1, p1, r1)
    y1, p1, r1 = 90, 0, 0
    R1 = geom.RfromYPR(y1, p1, r1)
    t1 = np.array([0,0,0])
    
    # camera 2
    w2, h2, fov2 = 640, 480, 49*np.pi/180
    K2 = geom.cameramat(w2,h2,fov2)
    #y2, p2, r2 = 100, 10, 3
    #R2 = geom.Rfrom_ypr(y2, p2, r2)
    #t2 = np.array([1,3,9])
    y2, p2, r2 = 150, 0, 0
    R2 = geom.RfromYPR(y2, p2, r2)
    t2 = np.array([8,0,4])
    
    # plane information and 3d points
#    n = np.array([-9,0,1])
#    offset = 87
#    pts3d = np.array([[10,0,3],
#                  [9.9,-1,2.1],
#                  [10.1,-0.5,3.9],
#                  [10.3,0.5,5.7]])
    n = np.array([1,0,-1])
    offset = -10
    inl = 20
    outl = 0
    inlrand = .1
    outrand = 4
    zrand = rand.normal(0,2,inl+outl)
    yrand = rand.normal(0,2,inl+outl)
    orand = np.array([rand.normal(0,outrand,3)+np.array([10,0,0]) for i in range(inl+outl)])
    orand[0:inl,:] = [(np.zeros(3) if inlrand==0 else rand.normal(0,inlrand,3)) for i in range(inl)]
    pts3d = [orand[i,:]+np.array([10+zrand[i],yrand[i],zrand[i]]) for i in range(inl+outl)]

#    pts3d = [np.array([10.1,-0.1,0.1]),
#             np.array([9,1,-1]),
#             np.array([11,2,1]),
#             np.array([12,-1,2])]

    # 3d points
    print '3d world points'
    print np.array(pts3d)

    # 3d image coordinate points
    print 'Camera 1 frame'
    pts3d1 = np.array([np.dot(np.transpose(R1),x-t1) for x in pts3d])
    print pts3d1
    print 'Camera 2 frame'
    pts3d2 = np.array([np.dot(np.transpose(R2),x-t2) for x in pts3d])
    print pts3d2

    # camera points
    Kpts1 = [np.dot(K1,np.dot(np.transpose(R1),x-t1)) for x in pts3d]
    Kpts1 = [p/p[2] for p in Kpts1]
    print 'Image coordinates for camera 1'
    print np.array(Kpts1)
    Kpts2 = [np.dot(K2,np.dot(np.transpose(R2),x-t2)) for x in pts3d]
    Kpts2 = [p/p[2] for p in Kpts2]
    print 'Image coordinates for camera 2'
    print np.array(Kpts2)

    # camere frame "3d points"
    pts1 = [alg.solve(K1,p) for p in Kpts1]
    pts1 = [p for p in pts1]
    print 'Inverse K camera 1 frame'
    print np.array(pts1)
    pts2 = [alg.solve(K2,p) for p in Kpts2]
    pts2 = [p for p in pts2]
    print 'Inverse K camera 2 frame'
    print np.array(pts2)

    matches = {}
    matches['qray'] = np.array(pts2)
    matches['dray'] = np.array(pts1)
    matches['w3d'] = np.array(pts3d)
    matches['nmat'] = len(pts3d)
    matches['numq'] = len(pts3d)
    matches['qidx'] = range(len(pts3d)+1)
    matches['hvrf'] = 0
    matches['herr'] = -1 * np.ones(len(pts3d))
#
#    matches = []
#    add_rand = 4
#    for (q,d,w) in zip(*[pts2,pts1,pts3d]):
#        numrand = rand.randint(add_rand)
#        dbs = [ {'d2d': np.array([int(d[1]),int(d[2])]), \
#                          'dray': d, \
#                          'dprm': [0,0], \
#                          'd3d': w, \
#                          'nnd': 0 } ]
#        for i in xrange(numrand):
#            worldrand = rand.normal(0,outrand,3)+np.array([10,0,0])
#            dbrand = np.dot(np.transpose(R1),worldrand-t1)
#            dbrand = dbrand / dbrand[2]
#            db = {'d2d': np.array([0,0]), \
#                  'dray': dbrand, \
#                  'dprm': [0,0], \
#                  'd3d': worldrand, \
#                  'nnd': 0 }
#            dbs.append(db)
#        match = {'q2d': [0,0], \
#                 'qray': q, \
#                 'qprm': [0,0], \
#                 'nmat': 1+numrand, \
#                 'hmat': -1, \
#                 'herr': -1, \
#                 'db': dbs }
#        matches.append(match)

#    print 'Matches lists of dictionaries'
#    print matches

    # test constrained homography pose estimation
    matches, t, nbear = computePose.constrainedHomography(matches,R1,R2,maxerr=.05,maxiter=1)

    print 'Number of inliers: %d' % matches['hvrf']
    print 'Query location: ' + str(t)
    print 'Plane normal bearing: %.1f degrees' % np.mod(180/np.pi*nbear,360)
    
#    # compute homography
#    H = cv.CreateMat(3, 3, cv.CV_64F)
#    cv.SetZero(H)
#    inliersH = cv.CreateMat(1, len(pts1), cv.CV_8U)
#    cv.SetZero(inliersH)
#    pts_q = cv.CreateMat(len(pts1), 1, cv.CV_64FC2)
#    pts_db = cv.CreateMat(len(pts1), 1, cv.CV_64FC2)
#    for i in range(len(pts1)):
#        cv.Set2D(pts_db, i, 0, cv.Scalar(*pts1[i,:2]))
#        cv.Set2D(pts_q, i, 0, cv.Scalar(*pts2[i,:2]))
#    cv.FindHomography(pts_db, pts_q, H, method=0)
#    inliersH = np.asarray(inliersH)[0]
#    H = np.asarray(H)
#    print 'H'
#    print H
#    print 'Homography transformation'
#    print np.dot(H,np.transpose(pts1)) / np.transpose(pts2)
#    dRH,tqH,nH = geom.decomposeH(H)
#    #xyzH = np.compress(inliersH,xyz,0)
#    #q2dH = np.compress(inliersH,arr(q2d),0)
#    #db2dH = np.compress(inliersH,arr(db2d),0)
#    #RqH = np.dot(Rd,dRH)
#
#    # this should be valid homography matrix as well
#    dR = np.dot(np.transpose(R2),R1)
#    dt = np.dot(np.transpose(R1),t2)
#    dt2 = np.dot(np.transpose(dR),dt)
#    n1 = np.dot(np.transpose(R1),n)
#    n1 = n1 / alg.norm(n1)
#    d = 10 / np.sqrt(2)
#    G = dR + (1/d)*np.outer(dt2,n1)
#    print 'Theory N'
#    print n1
#    print 'Theory T'
#    print dt2
#    print 'Theory 1/dTN'
#    print (1/d)*np.outer(dt2,n1)
#    print 'Theory R'
#    print np.transpose(dR)
#    print 'Theory H'
#    print G
#    print 'Theory homography transformation'
#    print np.dot(G,np.transpose(pts1)) / np.transpose(pts2)
#
#
#    pts1 = np.array(pts1)
#    pts2 = np.array(pts2)
#    # test out T from known R algorithm
#    dR = np.dot(np.transpose(R1),R2)
#    t, foo, mask, err = computePose.TfromR(dR, pts1, pts2)
#    print 'T from R'
#    print t
#    print 'Error in T from R'
#    print err
#    print 'Inliers in T from R'
#    print sum(mask)
