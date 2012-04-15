# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="aaronh"
__date__ ="$Nov 22, 2011 2:43:34 AM$"

import time
import scipy.optimize as opt
import numpy.random as rand
import numpy as np
import numpy.linalg as alg
from numpy import transpose as tp
import geom


def constrainedHomography(matches, wRd, wRq, qYaw=np.nan, nYaw=np.nan, runflag=0, maxerr=.05, iterscale=1000, yrestrict=False):

    print 'Solving constrained homography [%d]: %.3f starting error threshold, e-fold increase every %d iterations...' % (runflag,maxerr,iterscale)
    start = time.time()

#    print '*****************: %d, %d, %f' % (runflag,qYaw,nYaw)

    # Homography parameters to solve for:
    #   translation: 2 parameters (never known, always solved)
    #   normal yaw: 1 parameter (may be known. always solved)
    #   query yaw: 1 parameter (may be known. always solved)
    #   scale factor: 1 parameter (never known, may not solve for this)

    # Set the different run conditions to be used in the RANSAC loop
    # Note that if qYaw is unknown, then wRq is assumed just pitch and roll (yaw=0)
    if runflag == 0: # qYaw unknown, nYaw unknown
        nprm, nrand = 5, 3
        compH, lsqH, errH, confH = compH_tqn, lsqH_tqn, errH_tqn, confH_
    elif runflag == 1: # qYaw unknown, nYaw known
        nprm, nrand = 4, 2
        compH, lsqH, errH, confH = compH_tq, lsqH_tq, errH_tq, confH_
    elif runflag == 2: # qYaw known, nYaw unknown
        nprm, nrand = 4, 2
        compH, lsqH, errH, confH = compH_tn, lsqH_tn, errH_tn, confH_
    elif runflag == 3: # qYaw known, nYaw known
        nprm, nrand = 3, 2
        compH, lsqH, errH, confH = compH_t, lsqH_t, errH_t, confH_
    elif runflag == 4: # qYaw unknown, nYaw unknown, solve for depth
        nprm, nrand = 6, 3
        compH, lsqH, errH, confH = compH_dtqn, lsqH_dtqn, errH_dtqn, confH_d
    elif runflag == 5: # qYaw unknown, nYaw known, solve for depth
        nprm, nrand = 5, 2
        compH, lsqH, errH, confH = compH_dtq, lsqH_dtq, errH_dtq, confH_d
    elif runflag == 6: # qYaw known, nYaw unknown, solve for depth
        nprm, nrand = 5, 2
        compH, lsqH, errH, confH = compH_dtn, lsqH_dtn, errH_dtn, confH_d
    else: # runflag == 7: qYaw known, nYaw known, solve for depth
        nprm, nrand = 4, 2
        compH, lsqH, errH, confH = compH_dt, lsqH_dt, errH_dt, confH_d

    # Set variables
    nmat, numq = matches['nmat'], matches['numq']
    constants = (wRq,wRd,qYaw,nYaw)
    bprm, bmask, bnumi, bconf, berr, binlerr = np.zeros(nprm), np.bool_(np.zeros(nmat)), 0, 0, np.inf, np.nan
    qray, dray, ddep, qidx, weights = matches['qray'], matches['dray'], matches['ddep'], matches['qidx'], matches['weight']

    # Ransac loop to eliminate outliers with homography
    # Solves homography matrix for homography matrix H=qRd(I+tn') using y ~ Hx
    iter, maxiter, vcount, vstop, mininl = 0, 2*iterscale, 0, int(np.ceil(0.005*iterscale)), 15
    while vcount < vstop and iter < maxiter:
        iter += 1
        q, d, dep = randsamples(nrand, nmat, qray, dray, ddep)
        prm, valid = compH(q,d,dep,constants)
        if not valid: continue
        errs = errH(prm,qray,dray,ddep,constants)
        inlerr = maxerr*np.e**(float(iter)/iterscale)
        imask = getInliers(errs,inlerr,qidx,numq,nmat)
        ierr = rperr(errH,imask,prm,qray,dray,ddep,constants)
        numi, hconf = confH(prm,weights,imask,ierr,mininl,yrestrict)
        if hconf>0: vcount += 1
        if hconf>bconf: bprm, bmask, bnumi, bconf, berr, binlerr = prm, imask, numi, hconf, ierr, inlerr
    print 'Total iterations / valid solutions: %d | %d' % (iter,vcount)

    # Guided matching
    numi, imask, prm, hconf, ierr, inlerr = bnumi, bmask, bprm, bconf, berr, binlerr
    last_numi, gmiter, maxgm = 0, 0, 100
    while hconf > 0 and last_numi != numi and iter < maxgm:
        last_numi, gmiter = numi, gmiter+1
        q, d, dep = qray[imask,:], dray[imask,:], ddep[imask]
        prm = lsqH(prm,q,d,dep,constants)
        errs = errH(prm,qray,dray,ddep,constants)
        imask = getInliers(errs,inlerr,qidx,numq,nmat)
        ierr = rperr(errH,imask,prm,qray,dray,ddep,constants)
        numi, hconf = confH(prm,weights,imask,ierr,mininl,yrestrict)

    # Set output parameters
    matches['iter'] = iter
    matches['nvalid'] = vcount
    matches['inlerr'] = inlerr
    matches['iprm'] = prm
    matches['imask'] = imask
    matches['rperr'] = ierr
    matches['numi'] = numi
    matches['ifrc'] = float(numi) / matches['nmat']
    matches['hconf'] = 0 if matches['rperr']==0 else hconf / matches['rperr']
    matches['runflag'] = runflag

#    if matches['numi'] >= 15 and np.isnan(ierr):
#        print '---------------------------------------'
#        print imask
#        print np.sum(imask)
#        print ierr
#        print rperr(errH,imask,prm,qray,dray,ddep,constants)
#        print '---------------------------------------'

    # Print output state
    if vcount < vstop:
        print 'Constrained homography failed.'
        pose = np.zeros(6)
        pose[3:5] = np.nan
    else:
        print 'Resulting confidence of homography solution: %f' % matches['hconf']
        print 'Reprojection error / Error threshold of solution: %.3f / %.3f' % (matches['rperr'],matches['inlerr'])
        print 'Number of inliers / total correspondences: ' + str(matches['numi']) + ' / ' + str(nmat)
        if   np.mod(runflag,4) == 0: qYaw, nYaw = prm[3], prm[4]
        elif np.mod(runflag,4) == 1: qYaw, nYaw = prm[3], nYaw
        elif np.mod(runflag,4) == 2: qYaw, nYaw = qYaw, prm[3]
        if   runflag == 4: scale = prm[5]*geom.vecnorm(prm[:3])
        elif runflag == 5: scale = prm[4]*geom.vecnorm(prm[:3])
        elif runflag == 6: scale = prm[4]*geom.vecnorm(prm[:3])
        elif runflag == 7: scale = prm[3]*geom.vecnorm(prm[:3])
        else:              scale = np.nan
        pose = np.append( geom.normalrows(prm[:3]) , [ qYaw , nYaw , scale ] )
    print 'Constrained homography took %.1f seconds.' % (time.time()-start)

    return matches, pose

### global scale factor for plane distance error compared to reprojection error ###
### Equates 5% depth deviation to 1% lateral deviation ###
def alpha():
    return 0.2 # 5% depth error = 1% pixel reprojection error

def rperr(errH,imask,prm,qray,dray,ddep,constants):
    return geom.vecnorm(errH(prm,qray[imask,:],dray[imask,:],ddep[imask],constants)) / np.sum(imask)**0.5

def randsamples(nrand, nmat, qray, dray, ddep):
    i1 = rand.randint(0,nmat)
    q1, d1, dep1 = qray[i1,:], dray[i1,:], ddep[i1]
    i2 = rand.randint(0,nmat-1)
    i2 = i2+1 if i2>=i1 else i2
    q2, d2, dep2 = qray[i2,:], dray[i2,:], ddep[i2]
    if nrand == 2:
        q, d, dep = np.array([q1,q2]), np.array([d1,d2]), np.array([dep1,dep2])
    else: # nrand == 3
        i3 = rand.randint(0,nmat-2)
        i3 = i3 if i3<min(i1,i2) else ( i3+1 if i3+1<max(i1,i2) else i3+2 )
        q3, d3, dep3 = qray[i3,:], dray[i3,:], ddep[i3]
        q, d, dep = np.array([q1,q2,q3]), np.array([d1,d2,d3]), np.array([dep1,dep2,dep3])
    return q, d, dep

def getInliers(errs,maxerr,qidx,numq,nmat):
    emask = errs<maxerr
    imask = np.bool_(np.zeros(nmat))
    for i in xrange(numq):
        if not emask[qidx[i]:qidx[i+1]].any(): continue
        imask[qidx[i]+np.argmin(errs[qidx[i]:qidx[i+1]])] = True
    return imask

#####  RETURN CONFIDENCE OF HOMOGRAPHY ESTIMATE  #####

def confH_(prm,weights,imask,ierr,mininl,yrestrict):
    numi = np.sum(imask)
    iconf = np.sum(weights[imask]) / np.sum(weights)
    if numi<mininl: hconf = 0
    else: hconf = iconf #/ ierr
    return numi, hconf

def confH_d(prm,weights,imask,ierr,mininl,yrestrict):
    dist = prm[-1]*geom.vecnorm(prm[:3])
    ty = prm[-1]*prm[1]
    numi = np.sum(imask)
    ydev = abs(prm[-1]*prm[1]-2)
    yconf = 1 if ydev<2 or not yrestrict else 2/ydev
    iconf = np.sum(weights[imask]) / np.sum(weights)
    if numi<mininl or dist>50: hconf = 0
    else: hconf = yconf * iconf #/ ierr
    return numi, hconf

#####  COMPUTE HOMOGRAPHY MATRIX WITH MINIMUM NUMBER OF CORRESPONDENCES  #####

def compH_tqn(qray,dray,ddep,constants):
    # set variables
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    dRq = np.dot(tp(wRd),geom.RfromYPR(qYaw,pr[0],pr[1]))
    xd, yq = dray, qray
    yd = tp(np.dot(dRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yd,xd) # no renormalization to bias more confident planes
    # compute homography parameters
    teig = alg.eig(np.dot(tp(tn),tn))
    nullidx = np.argmin(teig[0])
    valid = teig[0][nullidx] < 1e-2
    t = geom.normalrows(teig[1][:,nullidx]) # homography translation
    m = geom.vecnorm(tn)/geom.vecnorm(np.cross(yd,t))*geom.vecnorm(xw[:,[0,2]])
    f = np.arctan2(xw[:,0],xw[:,2])
    errf = lambda prm,argm,argf: prm[0]-argm/np.cos(prm[1]-argf)
    kn_init = np.array([1.2*np.mean(m),np.mean(f)])
    k, n = tuple( opt.leastsq(errf,kn_init,args=(m,f),warning=False)[0] )
    valid = valid and np.std( m/(k*np.cos(n-f)) ) < 0.1
    fe = np.mod(n-np.mean(f),2*np.pi)
    if np.abs(fe) < np.pi/2: n = np.mod(n+np.pi,2*np.pi)
    if np.mean(np.inner(xd-yd,t)) < 0: t = -t
    # set parameters and refine
    prm = np.append(np.abs(k)*np.dot(wRd,t),[qYaw,180/np.pi*n])
    if valid: prm = lsqH_tqn(prm,qray,dray,ddep,constants)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_tq(qray,dray,ddep,constants):
    # set variables
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    dRq = np.dot(tp(wRd),geom.RfromYPR(qYaw,pr[0],pr[1]))
    xd, yq = dray, qray
    yd = tp(np.dot(dRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yd,xd)
    n = nYaw * np.pi/180 # homography normal bearing
    # compute homography parameters based off guessed yaw
    t = geom.normalrows(np.cross(tn[0,:],tn[1,:])) # homography translation
    m = geom.vecnorm(tn)/geom.vecnorm(np.cross(yd,t))*geom.vecnorm(xw[:,[0,2]])
    f = np.arctan2(xw[:,0],xw[:,2])
    k = np.mean( m / np.cos(n-f) )
    valid = np.std( m/(k*np.cos(n-f)) ) < 0.1
    fe = np.mod(n-np.mean(f),2*np.pi)
    if np.abs(fe) < np.pi/2: n = np.mod(n+np.pi,2*np.pi)
    if np.mean(np.inner(xd-yd,t)) < 0: t = -t
    # set parameters and refine
    prm = np.append(np.abs(k)*np.dot(wRd,t),qYaw)
    if valid: prm = lsqH_tq(prm,qray,dray,ddep,constants)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_tn(qray,dray,ddep,constants):
    # set variables
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    xd, yq = dray, qray
    yd = tp(np.dot(dRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yd,xd)
    # compute homography parameters
    t = geom.normalrows(np.cross(tn[0,:],tn[1,:])) # homography translation
    m = geom.vecnorm(tn)/geom.vecnorm(np.cross(yd,t))*geom.vecnorm(xw[:,[0,2]])
    f = np.arctan2(xw[:,0],xw[:,2])
    errf = lambda prm,argm,argf: prm[0]-argm/np.cos(prm[1]-argf)
    kn_init = np.array([1.2*np.mean(m),np.mean(f)])
    k, n = tuple( opt.leastsq(errf,kn_init,args=(m,f),warning=False)[0] )
    valid = np.std( m/(k*np.cos(n-f)) ) < 0.1
    fe = np.mod(n-np.mean(f),2*np.pi)
    if np.abs(fe) < np.pi/2: n = np.mod(n+np.pi,2*np.pi)
    if np.mean(np.inner(xd-yd,t)) < 0: t = -t
    # set parameters and refine
    prm = np.append(k*np.dot(wRd,t),180/np.pi*n)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_t(qray,dray,ddep,constants):
    # set variables
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    xd, yq = dray, qray
    yd = tp(np.dot(dRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yd,xd)
    n = nYaw * np.pi/180 # homography normal bearing
    # compute homography parameters
    t = geom.normalrows(np.cross(tn[0,:],tn[1,:])) # homography translation
    m = geom.vecnorm(tn)/geom.vecnorm(np.cross(yd,t))*geom.vecnorm(xw[:,[0,2]])
    f = np.arctan2(xw[:,0],xw[:,2])
    k = np.mean( m / np.cos(n-f) )
    valid = np.std( m/(k*np.cos(n-f)) ) < 0.1
    if np.mean(np.inner(xd-yd,t)) < 0: t = -t
    # set parameters and refine
    prm = np.abs(k)*np.dot(wRd,t)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_dtqn(qray,dray,ddep,constants):
    # set variables
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    dRq = np.dot(tp(wRd),geom.RfromYPR(qYaw,pr[0],pr[1]))
    xd, yq = dray, qray
    yd = tp(np.dot(dRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yd,xd) # no renormalization to bias more confident planes
    # compute homography parameters
    teig = alg.eig(np.dot(tp(tn),tn))
    nullidx = np.argmin(teig[0])
    valid = teig[0][nullidx] < 1e-2
    t = geom.normalrows(teig[1][:,nullidx]) # homography translation
    m = geom.vecnorm(tn)/geom.vecnorm(np.cross(yd,t))*geom.vecnorm(xw[:,[0,2]])
    f = np.arctan2(xw[:,0],xw[:,2])
    errf = lambda prm,argm,argf: prm[0]-argm/np.cos(prm[1]-argf)
    kn_init = np.array([1.2*np.mean(m),np.mean(f)])
    k, n = tuple( opt.leastsq(errf,kn_init,args=(m,f),warning=False)[0] )
    valid = valid and np.std( m/(k*np.cos(n-f)) ) < 0.1
    fe = np.mod(n-np.mean(f),2*np.pi)
    if np.abs(fe) < np.pi/2: n = np.mod(n+np.pi,2*np.pi)
    if np.mean(np.inner(xd-yd,t)) < 0: t = -t
    # compute plane depth
    nd = -np.dot(tp(wRd),[np.sin(n),0,np.cos(n)])
    dep = ddep*np.inner(dray,nd)
    pd = np.mean(dep)
    valid = valid and np.std(dep/pd) < 0.1
    # set parameters and refine
    prm = np.append(np.abs(k)*np.dot(wRd,t),[qYaw,180/np.pi*n,pd])
    if valid: prm = lsqH_dtqn(prm,qray,dray,ddep,constants)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_dtq(qray,dray,ddep,constants):
    # set variables
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    dRq = np.dot(tp(wRd),geom.RfromYPR(qYaw,pr[0],pr[1]))
    xd, yq = dray, qray
    yd = tp(np.dot(dRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yd,xd)
    n = nYaw * np.pi/180 # homography normal bearing
    # compute homography parameters based off guessed yaw
    t = geom.normalrows(np.cross(tn[0,:],tn[1,:])) # homography translation
    m = geom.vecnorm(tn)/geom.vecnorm(np.cross(yd,t))*geom.vecnorm(xw[:,[0,2]])
    f = np.arctan2(xw[:,0],xw[:,2])
    k = np.mean( m / np.cos(n-f) )
    valid = np.std( m/(k*np.cos(n-f)) ) < 0.1
    fe = np.mod(n-np.mean(f),2*np.pi)
    if np.abs(fe) < np.pi/2: n = np.mod(n+np.pi,2*np.pi)
    if np.mean(np.inner(xd-yd,t)) < 0: t = -t
    # compute plane depth
    nd = -np.dot(tp(wRd),[np.sin(n),0,np.cos(n)])
    dep = ddep*np.inner(dray,nd)
    pd = np.mean(dep)
    valid = valid and np.std(dep/pd) < 0.1
    # set parameters and refine
    prm = np.append(np.abs(k)*np.dot(wRd,t),[qYaw,pd])
    if valid: prm = lsqH_dtq(prm,qray,dray,ddep,constants)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_dtn(qray,dray,ddep,constants):
    # set variables
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    xd, yq = dray, qray
    yd = tp(np.dot(dRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yd,xd)
    # compute homography parameters
    t = geom.normalrows(np.cross(tn[0,:],tn[1,:])) # homography translation
    m = geom.vecnorm(tn)/geom.vecnorm(np.cross(yd,t))*geom.vecnorm(xw[:,[0,2]])
    f = np.arctan2(xw[:,0],xw[:,2])
    errf = lambda prm,argm,argf: prm[0]-argm/np.cos(prm[1]-argf)
    kn_init = np.array([1.2*np.mean(m),np.mean(f)])
    k, n = tuple( opt.leastsq(errf,kn_init,args=(m,f),warning=False)[0] )
    valid = np.std( m/(k*np.cos(n-f)) ) < 0.1
    fe = np.mod(n-np.mean(f),2*np.pi)
    if np.abs(fe) < np.pi/2: n = np.mod(n+np.pi,2*np.pi)
    if np.mean(np.inner(xd-yd,t)) < 0: t = -t
    # compute plane depth
    nd = -np.dot(tp(wRd),[np.sin(n),0,np.cos(n)])
    dep = ddep*np.inner(dray,nd)
    pd = np.mean(dep)
    valid = valid and np.std(dep/pd) < 0.1
    # set parameters and refine
    prm = np.append(k*np.dot(wRd,t),[180/np.pi*n,pd])
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_dt(qray,dray,ddep,constants):
    # set variables
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    xd, yq = dray, qray
    yd = tp(np.dot(dRq,tp(yq)))
    xw = tp(np.dot(wRd,tp(xd)))
    tn = np.cross(yd,xd)
    n = nYaw * np.pi/180 # homography normal bearing
    # compute homography parameters
    t = geom.normalrows(np.cross(tn[0,:],tn[1,:])) # homography translation
    m = geom.vecnorm(tn)/geom.vecnorm(np.cross(yd,t))*geom.vecnorm(xw[:,[0,2]])
    f = np.arctan2(xw[:,0],xw[:,2])
    k = np.mean( m / np.cos(n-f) )
    valid = np.std( m/(k*np.cos(n-f)) ) < 0.1
    if np.mean(np.inner(xd-yd,t)) < 0: t = -t
    # compute plane depth
    nd = -np.dot(tp(wRd),[np.sin(n),0,np.cos(n)])
    dep = ddep*np.inner(dray,nd)
    pd = np.mean(dep)
    valid = valid and np.std(dep/pd) < 0.1
    # set parameters and refine
    prm = np.append(k*np.dot(wRd,t),pd)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

##### SOLVE FOR HOMOGRAPHY MATRIX THAT MINIMIZES REPROJECTION ERROR #####

def lsqH_tqn(prm,qray,dray,ddep,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    try: return opt.leastsq(homerrf_tqn,prm,args=(qray,dray,ddep,pr,wRd),warning=False)[0]
    except TypeError: return prm

def lsqH_tq(prm,qray,dray,ddep,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return opt.leastsq(homerrf_tq,prm,args=(qray,dray,ddep,pr,wRd,nYaw),warning=False)[0]

def lsqH_tn(prm,qray,dray,ddep,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return opt.leastsq(homerrf_tn,prm,args=(qray,dray,ddep,dRq,wRd),warning=False)[0]

def lsqH_t(prm,qray,dray,ddep,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return opt.leastsq(homerrf_t,prm,args=(qray,dray,ddep,dRq,wRd,nYaw),warning=False)[0]

def lsqH_dtqn(prm,qray,dray,ddep,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    try: return opt.leastsq(homerrf_dtqn,prm,args=(qray,dray,ddep,pr,wRd),warning=False)[0]
    except TypeError: return prm

def lsqH_dtq(prm,qray,dray,ddep,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return opt.leastsq(homerrf_dtq,prm,args=(qray,dray,ddep,pr,wRd,nYaw),warning=False)[0]

def lsqH_dtn(prm,qray,dray,ddep,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return opt.leastsq(homerrf_dtn,prm,args=(qray,dray,ddep,dRq,wRd),warning=False)[0]

def lsqH_dt(prm,qray,dray,ddep,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return opt.leastsq(homerrf_dt,prm,args=(qray,dray,ddep,dRq,wRd,nYaw),warning=False)[0]

#####  COMPUTE CORRESPONDENCE ERRORS FOR HOMOGRAPHY COMPUTATIONS  #####

def homerrf_tqn(prm,qray,dray,ddep,pr,wRd):
    # set variables
    dRq = np.dot(tp(wRd),geom.RfromYPR(prm[3],pr[0],pr[1]))
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(prm[4]*np.pi/180),0,np.cos(prm[4]*np.pi/180)])
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.append( qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] , qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] )
    return err

def homerrf_tq(prm,qray,dray,ddep,pr,wRd,nbear):
    # set variables
    dRq = np.dot(tp(wRd),geom.RfromYPR(prm[3],pr[0],pr[1]))
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(nbear*np.pi/180),0,np.cos(nbear*np.pi/180)])
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.append( qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] , qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] )
    return err

def homerrf_tn(prm,qray,dray,ddep,dRq,wRd):
    # set variables
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(prm[3]*np.pi/180),0,np.cos(prm[3]*np.pi/180)])
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.append( qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] , qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] )
    return err

def homerrf_t(prm,qray,dray,ddep,dRq,wRd,nbear):
    # set variables
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(nbear*np.pi/180),0,np.cos(nbear*np.pi/180)])
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.append( qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] , qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] )
    return err

def homerrf_dtqn(prm,qray,dray,ddep,pr,wRd):
    # set variables
    dRq = np.dot(tp(wRd),geom.RfromYPR(prm[3],pr[0],pr[1]))
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(prm[4]*np.pi/180),0,np.cos(prm[4]*np.pi/180)])
    pd = prm[5]
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.concatenate( [ qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] ,
                            qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] ,
                            alpha() * ( pd - ddep*np.inner(dray,nd) ) / pd ] )
    return err

def homerrf_dtq(prm,qray,dray,ddep,pr,wRd,nbear):
    # set variables
    dRq = np.dot(tp(wRd),geom.RfromYPR(prm[3],pr[0],pr[1]))
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(nbear*np.pi/180),0,np.cos(nbear*np.pi/180)])
    pd = prm[4]
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.concatenate( [ qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] ,
                            qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] ,
                            alpha() * ( pd - ddep*np.inner(dray,nd) ) / pd ] )
    return err

def homerrf_dtn(prm,qray,dray,ddep,dRq,wRd):
    # set variables
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(prm[3]*np.pi/180),0,np.cos(prm[3]*np.pi/180)])
    pd = prm[4]
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.concatenate( [ qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] ,
                            qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] ,
                            alpha() * ( pd - ddep*np.inner(dray,nd) ) / pd ] )
    return err

def homerrf_dt(prm,qray,dray,ddep,dRq,wRd,nbear):
    # set variables
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(nbear*np.pi/180),0,np.cos(nbear*np.pi/180)])
    pd = prm[3]
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.concatenate( [ qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] ,
                            qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] ,
                            alpha() * ( pd - ddep*np.inner(dray,nd) ) / pd ] )
    return err

#####  COMPUTE OVERALL ERROR FOR CORRESPONDENCES  #####

def errH_tqn(prm,qray,dray,ddep,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return np.sqrt(np.sum(np.reshape(homerrf_tqn(prm,qray,dray,ddep,pr,wRd),[2,-1])**2,0))

def errH_tq(prm,qray,dray,ddep,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return np.sqrt(np.sum(np.reshape(homerrf_tq(prm,qray,dray,ddep,pr,wRd,nYaw),[2,-1])**2,0))

def errH_tn(prm,qray,dray,ddep,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return np.sqrt(np.sum(np.reshape(homerrf_tn(prm,qray,dray,ddep,dRq,wRd),[2,-1])**2,0))

def errH_t(prm,qray,dray,ddep,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return np.sqrt(np.sum(np.reshape(homerrf_t(prm,qray,dray,ddep,dRq,wRd,nYaw),[2,-1])**2,0))

def errH_dtqn(prm,qray,dray,ddep,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return np.sqrt(np.sum(np.reshape(homerrf_dtqn(prm,qray,dray,ddep,pr,wRd),[3,-1])**2,0))

def errH_dtq(prm,qray,dray,ddep,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return np.sqrt(np.sum(np.reshape(homerrf_dtq(prm,qray,dray,ddep,pr,wRd,nYaw),[3,-1])**2,0))

def errH_dtn(prm,qray,dray,ddep,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return np.sqrt(np.sum(np.reshape(homerrf_dtn(prm,qray,dray,ddep,dRq,wRd),[3,-1])**2,0))

def errH_dt(prm,qray,dray,ddep,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return np.sqrt(np.sum(np.reshape(homerrf_dt(prm,qray,dray,ddep,dRq,wRd,nYaw),[3,-1])**2,0))