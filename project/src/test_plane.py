# To change this template, choose Tools | Templates
# and open the template in the editor.

import numpy as np
import numpy.random as rnd
import numpy.linalg as alg
import geom
import computePose

if __name__ == '__main__':

    # plane prms
    n = geom.normalrows(rnd.randn(3))
    n[1] *= 0.1
    d = 10*np.abs(rnd.randn())
    print 'Normal vector: [ %.2f , %.2f , %.2f ]' % (n[0],n[1],n[2])
    print 'Plane offset: %.1f' % d
    
    # 3d pts
    npts = 2000
    ninl = 500
    x = 10*rnd.randn(npts,1)
    y = 10*rnd.randn(npts,1)
    z = np.zeros([npts,1])
#    print np.shape(z[:ninl,0])
#    print np.shape(rnd.randn(ninl,1))
#    print np.shape(x[:ninl,0])
    z[:ninl,0] = 1/n[2] * ( d + 0*rnd.randn(ninl) - n[0]*x[:ninl,0] - n[1]*y[:ninl,0] )
    z[ninl:,0] = 1/n[2] * ( d + 20*rnd.randn(npts-ninl) - n[0]*x[ninl:,0] - n[1]*y[ninl:,0] )
    #z = 1/n[2] * (d-n[0]*x-n[1]*y)
    pts = np.concatenate( (x,y,z,np.ones([npts,1])) , 1 )
    xz_pts = pts[:,[0,2,3]]

    # RANSAC solve
    threshold = 1 # meters
    bprm, bnumi, bconf, bmask = np.zeros(4), 0, 0, np.bool_(np.zeros(npts))
    for i in range(10000):
        i1 = rnd.randint(0,npts)
        i2 = rnd.randint(0,npts-1)
        i2 = i2 if i2<i1 else i2+1
        i3 = rnd.randint(0,npts-2)
        i3 = i3 if i3<min(i1,i2) else ( i3+1 if i3+1<max(i1,i2) else i3+2 )
        inlpts = xz_pts[[i1,i2,i3],:]
        prm = geom.smallestSingVector(inlpts)
        prm = prm / geom.vecnorm(prm[:2])
        prm = prm if prm[2]<0 else -prm
        errs = np.abs(np.inner(xz_pts,prm))
        inlmask = errs < threshold
        numi = np.sum(inlmask)
        conf = numi**1.5 / np.sum(errs[inlmask]) * ( numi > 20 )
        if conf > bconf: bprm, bmask, bconf, bnumi = prm, inlmask, conf, numi
    prm, numi, mask = bprm, bnumi, bmask
    print '--- AFTER RANSAC ---'
    print 'Number of inliers / number of points : %d / %d' % (numi,npts)
    print 'Computed normal vector: [ %.2f , %.2f , %.2f ]' % (prm[0],0,prm[1])
    print 'Computed plane offset: %.1f' % -prm[2]

    # RANSAC solve
    threshold = 1 # meters
    bprm, bnumi, bconf, bmask = np.zeros(4), 0, 0, np.bool_(np.zeros(npts))
    for i in range(10000):
        i1 = rnd.randint(0,npts)
        i2 = rnd.randint(0,npts-1)
        i2 = i2 if i2<i1 else i2+1
        i3 = rnd.randint(0,npts-2)
        i3 = i3 if i3<min(i1,i2) else ( i3+1 if i3+1<max(i1,i2) else i3+2 )
        inlpts = xz_pts[[i1,i2,i3],:]
        prm = geom.smallestSingVector(inlpts)
        prm = prm / geom.vecnorm(prm[:2])
        prm = prm if prm[2]<0 else -prm
        errs = np.abs(np.inner(xz_pts,prm))
        inlmask = errs < threshold
        numi = np.sum(inlmask)
        conf = numi
        if conf > bconf: bprm, bmask, bconf, bnumi = prm, inlmask, conf, numi
    prm, numi, mask = bprm, bnumi, bmask
    print '--- AFTER RANSAC ---'
    print 'Number of inliers / number of points : %d / %d' % (numi,npts)
    print 'Computed normal vector: [ %.2f , %.2f , %.2f ]' % (prm[0],0,prm[1])
    print 'Computed plane offset: %.1f' % -prm[2]

    # guided matching
    for i in range(10):
        prm = geom.smallestSingVector(xz_pts[mask,:])
        prm = prm / geom.vecnorm(prm[:2])
        prm = prm if prm[2]<0 else -prm
        errs = np.abs(np.inner(xz_pts,prm))
        mask = errs < threshold
        numi = np.sum(mask)
    print '--- AFTER GUIDED MATCHING ---'
    print 'Number of inliers / number of points : %d / %d' % (numi,npts)
    print 'Computed normal vector: [ %.2f , %.2f , %.2f ]' % (prm[0],0,prm[1])
    print 'Computed plane offset: %.1f' % -prm[2]

    # Compute error
    print 'Plane deviation = %.2f meters' % np.mean(np.abs(np.inner(xz_pts[mask,:],prm)))