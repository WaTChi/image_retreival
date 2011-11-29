# test homography matrix solver

import solveHomography
import computePose
import numpy as np
import numpy.linalg as alg
import numpy.random as rnd
from numpy import transpose as tp
import geom
import computePose2

if __name__ == '__main__':

    # set useful constants and parameters
    numq = 500
    nperq = 3
    frac_inl = 0.3
    nout = int((1-frac_inl)*numq)
    ninl = numq-nout
    npts = numq*nperq
    zh = np.array([0,0,1])
    pi = np.pi

    # set up database parameters
    dyaw, dpitch, droll = 360*rnd.rand(), 10*rnd.rand()-5, 0
    #dyaw, dpitch, droll = 90, 0, 0
    dwidth, dheight, dfov = 2500, 1200, 90*pi/180
    Kd = geom.cameramat(dwidth,dheight,dfov)
    Kdinv = alg.inv(Kd)
    wRd = geom.RfromYPR(dyaw,dpitch,droll)
    dbcenter = np.dot(wRd,zh)

    # set up scene parameters
    depth = 10+5*rnd.randn()**2
    #depth = 10
    center = depth*dbcenter
    pyaw, ppitch, proll = np.mod( dyaw+180+60*rnd.rand()-30 , 360 ) , 0 , 0
    #pyaw = -45
    pnorm = np.dot( geom.RfromYPR(pyaw+180,ppitch,proll) , zh )
    plane_depth = np.inner(pnorm,center)
    peq = np.append( pnorm , -plane_depth )

    # set up query parameters
    qdepth = 10+5*rnd.randn()**2
    #qdepth = 10
    qcenter = center + 0.2*depth*(rnd.rand(3)-0.5)
    #qcenter = center
    qyaw, qpitch, qroll = np.mod( dyaw+60*rnd.rand()-30 , 360 ) , 20*rnd.rand()-10 , 0
    #qyaw, qpitch, qroll = 180, 0, 0
    qwidth, qheight, qfov = 1600, 1200, 60*pi/180
    Kq = geom.cameramat(qwidth,qheight,qfov)
    Kqinv = alg.inv(Kq)
    wRq = geom.RfromYPR(qyaw,qpitch,qroll)
    qloc = qcenter - qdepth*np.dot(wRq,zh)

    # generate feature points
    ngen = numq*nperq+nout
    dray = geom.normalrows( tp( np.dot( Kdinv , tp( np.concatenate( ( 500+1500*rnd.rand(ngen,1) , 300+600*rnd.rand(ngen,1), np.ones((ngen,1)) ) , 1 ) ) ) ) )
    w3d = tp( np.dot( wRd , tp(dray) ) )
    scale = plane_depth / np.inner(w3d,pnorm)
    w3d = tp( tp(w3d) * np.tile(scale,[3,1]) )
    qray = geom.normalrows( tp( np.dot( tp(wRq) , tp( w3d - np.tile(qloc,[ngen,1]) ) ) ) )
    for i in range(numq):
        if i<ninl:
            qray[range(nperq*i,nperq*(i+1)),:] = np.tile( qray[nperq*i,:] , [nperq,1] )
        else:
            qray[range(nperq*i,nperq*(i+1)),:] = np.tile( qray[npts+numq-i-1,:] , [nperq,1] )
    qray, dray, w3d = qray[:npts,:], dray[:npts,:], w3d[:npts,:]
    
    # print scene geometry
    print 'Query location : %.1f , %.1f, %.1f' % (qloc[0],qloc[1],qloc[2])
    print 'Query YPR      : %.1f , %.1f, %.1f' % (qyaw,qpitch,qroll)
    print 'Database YPR   : %.1f , %.1f, %.1f' % (dyaw,dpitch,droll)
    print 'Plane yaw      : %.1f' % pyaw

#    # print correspondences
#    print '3d points:'
#    print w3d
#    print 'Database rays:'
#    print dray
#    print 'Query rays:'
#    print qray
    
    # check to see if points are valid
    inlidx = range(0,nperq*ninl,nperq)
    print inlidx
    qvalid, dvalid, w3dvalid = qray[inlidx,:], dray[inlidx,:], w3d[inlidx,:]
    prm = qloc / plane_depth
    prm_tn = np.append( prm , pyaw )
    prm_tq = np.append( prm , qyaw )
    prm_tqn = np.append( prm , [qyaw,pyaw] )
    constants = (wRq,wRd,qyaw,pyaw)
    print 'Error of theory Homography :'
    print np.mean(solveHomography.errH_t(prm,qvalid,dvalid,constants))
    print np.mean(solveHomography.errH_tn(prm_tn,qvalid,dvalid,constants))
    print np.mean(solveHomography.errH_tq(prm_tq,qvalid,dvalid,constants))
    print np.mean(solveHomography.errH_tqn(prm_tqn,qvalid,dvalid,constants))
    prm2 = np.append( prm , [ np.pi + pyaw*np.pi/180 , plane_depth ] )
    dRq = np.dot( tp(wRd) , wRq )
    print 1.5*np.mean(computePose2.homerrf_RNP(prm2, qvalid, dvalid, w3d, dRq, wRd, 0))

    # set up matches dictionary
    matches = {}
    matches['qray'] = qray
    matches['dray'] = dray
    matches['w3d'] = w3d
    matches['qidx'] = np.arange(0,nperq*(numq+1),nperq)
    matches['nmat'] = npts
    matches['numq'] = numq

    # solve homography and extract parameters and compute scale
    matches, pose = solveHomography.constrainedHomography(matches,wRd,wRq,qyaw,pyaw,0,.03,100)
    tray, qbear, pbear = pose[:3], pose[3], pose[4]

    matches2, tray2, foo = computePose2.constrainedHomography(matches,wRd,wRq,np.nan,.03,100)

    # compute scaled translation from 3d pts
    comp_wRq = geom.RfromYPR(qbear,qpitch,qroll)
    t, scale, sdev = computePose.scalefrom3d(matches, tray, comp_wRq)

    t2, scale2, sdev2 = computePose.scalefrom3d(matches2, tray2, wRq)

    # print computed pose
    print 'Computed query location : %.1f , %.1f, %.1f' % (t[0],t[1],t[2])
    print 'Computed query yaw      : %.1f' % qbear
    print 'Computed plane yaw      : %.1f' % pbear
    print 'Mean scale              : %.1f' % scale
    print 'Scale std dev           : %.1f' % sdev
    
    # print computer pose 2
    print '[2]Computed query location : %.1f , %.1f, %.1f' % (t2[0],t2[1],t2[2])
    print '[2]Mean scale              : %.1f' % scale2
    print '[2]Scale std dev           : %.1f' % sdev2



    print 'Done.'

