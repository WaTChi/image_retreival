# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="aaronh"
__date__ ="$Nov 22, 2011 2:43:34 AM$"

import time
import scipy.optimize.minpack as opt
import numpy.random as rand
import numpy as np
import numpy.linalg as alg
from numpy import transpose as tp
import geom


def constrainedEssMatrix(matches, wRd, wRq, qYaw=np.nan, runflag=0, maxerr=.05, maxiter=1000):

    print 'Solving constrained essential matrix...'
    start = time.time()

    # Homography parameters to solve for:
    #   translation: 2 parameters (never known)
    #   normal yaw: 1 parameter (may be known)
    #   query yaw: 1 parameter (may be known)

    # Set the different run conditions to be used in the RANSAC loop
    # Note that if qYaw is unknown, then wRq is assumed just pitch and roll (yaw=0)
    if runflag == 0: # qYaw unknown
        nprm, nrand = 3, 3
        compE, lsqE, errE = compE_tq, lsqE_tq, errE_tq
    else: # runflag == 1: qYaw known
        nprm, nrand = 2, 2
        compE, lsqE, errE = compE_t, lsqE_t, errE_t

    # Set variables
    nmat, numq = matches['nmat'], matches['numq']
    constants = (wRq,wRd,qYaw)
    bprm, bmask, bnumi, bdomidx, berr = np.zeros(nprm), np.zeros(nmat), 0, -1, np.inf
    qray, dray, qidx = matches['qray'], matches['dray'], matches['qidx']

    # Ransac loop to eliminate outliers with essential matrix
    # Solves essential matrix
    iter = 0
    while iter < maxiter:
        iter += 1
        q, d = randsamples(nrand, nmat, qray, dray)
        prm, domidx, valid = compE(q,d,constants)
        if not valid: continue
        errs = errE(prm,qray,dray,constants,domidx)
        imask, numi = getInliers(errs,maxerr,qidx,numq,nmat)
        if numi >= bnumi: bprm, bmask, bnumi, bdomidx = prm, imask, numi, domidx
    # Guided matching
    numi, imask, prm, domidx = bnumi, bmask, bprm, bdomidx
    last_numi, iter, maxgm = 0, 0, 100
    while last_numi != numi and iter < maxgm:
        last_numi, iter = numi, iter+1
        q, d = qray[imask,:], dray[imask,:]
        prm = lsqE(prm,q,d,constants,domidx)
        errs = errE(prm,qray,dray,constants,domidx)
        imask, numi = getInliers(errs,maxerr,qidx,numq,nmat)

    # Set output parameters
    matches['iprm'] = prm
    matches['imask'] = imask
    matches['ierr'] = geom.vecnorm(errE(prm,qray[imask,:],dray[imask,:],constants,domidx)) * numq / numi if numi!=0 else np.inf
    matches['numi'] = sum(imask)

    # Print output state
    if matches['numi'] == 0:
        print 'Constrained homography failed.'
        pose = np.nan * np.zeros(4)
    else:
        print 'Result from error metric choosing best inlier set: %f' % matches['ierr']
        print 'Number of inliers / total correspondences: ' + str(matches['numi']) + ' / ' + str(nmat)
        if runflag == 0: qYaw = prm[3]
        pose = np.append( geom.normalrows(prm[:3]) , qYaw )
    print 'Constrained homography took %.1f seconds.' % (time.time()-start)

    return matches, pose


def randsamples(nrand, nmat, qray, dray):
    i1 = rand.randint(0,nmat)
    q1, d1 = qray[i1,:], dray[i1,:]
    i2 = rand.randint(0,nmat-1)
    i2 = i2+1 if i2>=i1 else i2
    q2, d2 = qray[i2,:], dray[i2,:]
    if nrand == 2:
        q, d = np.array([q1,q2]), np.array([d1,d2])
    else: # nrand == 3
        i3 = rand.randint(0,nmat-2)
        i3 = ( i3+2 if i3>=max(i1,i2) else (i3+1 if i3>=min(i1,i2) else i3) )
        q3, d3 = qray[i3,:], dray[i3,:]
        q, d = np.array([q1,q2,q3]), np.array([d1,d2,d3])
    return q, d

def getInliers(errs,maxerr,qidx,numq,nmat):
    emask = errs<maxerr
    imask = np.bool_(np.zeros(nmat))
    for i in xrange(numq):
        if not emask[qidx[i]:qidx[i+1]].any(): continue
        imask[qidx[i]+np.argmin(errs[qidx[i]:qidx[i+1]])] = True
    return imask, np.sum(imask)

#####  COMPUTE ESSENTIAL MATRIX WITH MINIMUM NUMBER OF CORRESPONDENCES  #####

def compE_tq(qray,dray,constants):
    # set variables
    Rpr, wRd, qYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    wRq = geom.RfromYPR(qYaw,pr[0],pr[1])
    xd, yq = dray, qray
    yw = tp(np.dot(wRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yw,xw) # no renormalization to bias more confident planes
    # compute essential matrix parameters based off guessed yaw
    teig = alg.eig(np.dot(tp(tn),tn))
    nullidx = np.argmin(teig[0])
    valid = teig[0][nullidx]/teig[0][np.argmax(teig[0])] < 1e-2
    t = geom.normalrows(teig[1][:,nullidx]) # essential matrix translation
    domidx = np.argmax(t)
    prm = np.append( np.delete(t/t[domidx],domidx) , qYaw )
    if valid: prm = lsqE_tq(prm,qray,dray,constants,domidx)
    return prm, domidx, valid

def compE_t(qray,dray,constants):
    # set variables
    wRq, wRd, qYaw = constants
    xd, yq = dray, qray
    yw = tp(np.dot(wRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yw,xw)
    # compute essential matrix parameters based off guessed yaw
    t = geom.normalrows(np.cross(tn[0,:],tn[1,:])) # homography translation
    return t, -1, True

##### SOLVE FOR ESSENTIAL MATRIX THAT MINIMIZES REPROJECTION ERROR #####

def lsqE_tq(prm,qray,dray,constants,domidx):
    Rpr, wRd, qYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return opt.leastsq(esserrf_tq,prm,args=(qray,dray,pr,wRd,domidx),warning=False)[0]

def lsqE_t(prm,qray,dray,constants,domidx):
    # set variables
    wRq, wRd, qYaw = constants
    xd, yq = dray, qray
    yw = tp(np.dot(wRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yw,xw) # no renormalization to bias more confident planes
    # compute essential matrix parameters based off guessed yaw
    teig = alg.eig(np.dot(tp(tn),tn))
    return geom.normalrows(teig[1][:,np.argmin(teig[0])]) # essential matrix translation

#####  COMPUTE CORRESPONDENCE ERRORS FOR ESSENTIAL MAT COMPUTATIONS  #####

def esserrf_tq(prm,qray,dray,pr,wRd,domidx):
    # set variables
    dRq = np.dot(tp(wRd),geom.RfromYPR(prm[2],pr[0],pr[1]))
    td = np.dot(tp(wRd),geom.normalrows(np.insert(prm[:2],domidx,1)))
    E = np.dot(tp(dRq),geom.xprodmat(td))
    # Compute homography error
    return np.sum( qray * tp(np.dot(E,tp(dray))) , 1 )

#####  COMPUTE OVERALL ERROR FOR CORRESPONDENCES  #####

def errE_tq(prm,qray,dray,constants,domidx):
    Rpr, wRd, qYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return np.abs( esserrf_tq(prm,qray,dray,pr,wRd,domidx) )

def errE_t(prm,qray,dray,constants,domidx):
    # set variables
    wRq, wRd, qYaw = constants
    dRq = np.dot(tp(wRd),wRd)
    td = np.dot(tp(wRd),prm)
    E = np.dot(tp(dRq),geom.xprodmat(td))
    # Compute homography error
    return np.abs( np.sum( qray * tp(np.dot(E,tp(dray))) , 1 ) )