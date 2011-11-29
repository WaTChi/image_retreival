# test homography matrix solver

import solveEssMatrix
import numpy as np
import numpy.linalg as alg
import numpy.random as rnd
from numpy import transpose as tp
import geom

if __name__ == '__main__':

    # set useful constants and parameters
    npts = 3
    zh = np.array([0,0,1])
    pi = np.pi

    # set up database parameters
    dyaw, dpitch, droll = 360*rnd.rand(), 10*rnd.rand()-5, 0
    dyaw, dpitch, droll = 90, 0, 0
    dwidth, dheight, dfov = 2500, 1200, 90*pi/180
    Kd = geom.cameramat(dwidth,dheight,dfov)
    Kdinv = alg.inv(Kd)
    wRd = geom.RfromYPR(dyaw,dpitch,droll)
    dbcenter = np.dot(wRd,zh)

    # set up scene parameters
    depth = 10*rnd.randn()**2
    depth = 10
    center = depth*dbcenter
    pyaw, ppitch, proll = np.mod( dyaw+180+60*rnd.rand()-30 , 360 ) , 0 , 0
    pyaw = -45
    pnorm = np.dot( geom.RfromYPR(pyaw+180,ppitch,proll) , zh )
    plane_depth = np.inner(pnorm,center)
    peq = np.append( pnorm , -plane_depth )

    # set up query parameters
    qdepth = 10*rnd.randn()**2
    qdepth = 10
    qcenter = center + 0.2*depth*(rnd.rand(3)-0.5)
    qcenter = center
    qyaw, qpitch, qroll = np.mod( dyaw+60*rnd.rand()-30 , 360 ) , 20*rnd.rand()-10 , 0
    qyaw, qpitch, qroll = 180, 0, 0
    qwidth, qheight, qfov = 1600, 1200, 60*pi/180
    Kq = geom.cameramat(qwidth,qheight,qfov)
    Kqinv = alg.inv(Kq)
    wRq = geom.RfromYPR(qyaw,qpitch,qroll)
    qloc = qcenter - qdepth*np.dot(wRq,zh)

    # generate feature points
    dray = geom.normalrows( tp( np.dot( Kdinv , tp( np.concatenate( ( 1100+300*rnd.rand(npts,1) , 500+200*rnd.rand(npts,1), np.ones((npts,1)) ) , 1 ) ) ) ) )
    print dray
    w3d = tp( np.dot( wRd , tp(dray) ) )
    scale = plane_depth / np.inner(w3d,pnorm)
    print plane_depth
    print np.inner(w3d,pnorm)
    print scale
    w3d = tp( tp(w3d) * np.tile(scale,[3,1]) )
    qray = geom.normalrows( tp( np.dot( tp(wRq) , tp( w3d - np.tile(qloc,[npts,1]) ) ) ) )

    # print scene geometry
    print 'Query location : %.1f , %.1f, %.1f' % (qloc[0],qloc[1],qloc[2])
    print 'Query YPR      : %.1f , %.1f, %.1f' % (qyaw,qpitch,qroll)
    print 'Database YPR   : %.1f , %.1f, %.1f' % (dyaw,dpitch,droll)
    print 'Plane yaw      : %.1f' % pyaw
    print '3d points:'
    print w3d

    # check to see if points are valid
    prm = geom.normalrows(qloc)
    constants = (wRq,wRd,qyaw)
    print prm
    print 'Error of theory Essential matrix :'
    print solveEssMatrix.errE_t(prm,qray,dray,constants,-1)

    # set up matches dictionary
    matches = {}
    matches['qray'] = qray
    matches['dray'] = dray
    matches['w3d'] = w3d
    matches['qidx'] = np.arange(npts+1)
    matches['nmat'] = npts
    matches['numq'] = npts

    # solve homography
    matches = solveEssMatrix.constrainedEssMatrix(matches,wRd,wRq,qyaw,6,.01,1)

    print 'Done.'