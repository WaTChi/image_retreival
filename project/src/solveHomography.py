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


def constrainedHomography(matches, wRd, wRq, qYaw=np.nan, nYaw=np.nan, runflag=0, maxerr=.05, maxiter=1000):

    print 'Solving constrained homography...'
    start = time.time()

    # Homography parameters to solve for:
    #   translation: 3 parameters (never known)
    #   normal yaw: 1 parameter (may be known)
    #   query yaw: 1 parameter (may be known)

    # Set the different run conditions to be used in the RANSAC loop
    # Note that if qYaw is unknown, then wRq is assumed just pitch and roll (yaw=0)
    if runflag == 0: # qYaw unknown, nYaw unknown
        nprm, nrand = 5, 3
        compH, lsqH, errH = compH_tqn, lsqH_tqn, errH_tqn
    elif runflag == 1: # qYaw unknown, nYaw known
        nprm, nrand = 4, 2
        compH, lsqH, errH = compH_tq, lsqH_tq, errH_tq
    elif runflag == 2: # qYaw known, nYaw unknown
        nprm, nrand = 4, 2
        compH, lsqH, errH = compH_tn, lsqH_tn, errH_tn
    else: # runflag == 3: qYaw known, nYaw known
        nprm, nrand = 3, 2
        compH, lsqH, errH = compH_t, lsqH_t, errH_t

    # Set variables
    nmat, numq = matches['nmat'], matches['numq']
    constants = (wRq,wRd,qYaw,nYaw)
    bprm, bmask, bnumi, berr = np.zeros(nprm), np.bool_(np.zeros(nmat)), 0, np.inf
    qray, dray, qidx = matches['qray'], matches['dray'], matches['qidx']

    # Ransac loop to eliminate outliers with homography
    # Solves homography matrix for homography matrix H=qRd(I+tn') using y ~ Hx
    iter = 0
    while iter < maxiter:
        iter += 1
        q, d = randsamples(nrand, nmat, qray, dray)
        prm, valid = compH(q,d,constants)
        if not valid: continue
        errs = errH(prm,qray,dray,constants)
        imask, numi = getInliers(errs,maxerr,qidx,numq,nmat)
#        ierr = geom.vecnorm(errH(prm,qray[imask,:],dray[imask,:],constants)) * numq / numi
#        if ierr < berr: berr, bprm, bmask, bnumi = ierr, prm, imask, numi
        if numi >= bnumi: bprm, bmask, bnumi = prm, imask, numi

    # Guided matching
    numi, imask, prm = bnumi, bmask, bprm
    last_numi, iter, maxgm = 0, 0, 100
    while last_numi != numi and iter < maxgm:
        last_numi, iter = numi, iter+1
        q, d = qray[imask,:], dray[imask,:]
        prm = lsqH(prm,q,d,constants)
        errs = errH(prm,qray,dray,constants)
        imask, numi = getInliers(errs,maxerr,qidx,numq,nmat)

    # Set output parameters
    matches['iprm'] = prm
    matches['imask'] = imask
    matches['ierr'] = geom.vecnorm(errH(prm,qray[imask,:],dray[imask,:],constants)) * numq / numi
    matches['nvrf'] = sum(imask)

    # Print output state
    if matches['nvrf'] == 0:
        print 'Constrained homography failed.'
        pose = np.nan*np.zeros(5)
    else:
        print 'Result from error metric choosing best inlier set: %f' % matches['ierr']
        print 'Number of inliers / total correspondences: ' + str(matches['nvrf']) + ' / ' + str(nmat)
        if   runflag == 0: qYaw, nYaw = prm[3], prm[4]
        elif runflag == 1: qYaw = prm[3]
        elif runflag == 2: nYaw = prm[3]
        pose = np.append( geom.normalrows(prm[:3]) , [ qYaw , nYaw ] )
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
        i3 = i3 if i3<min(i1,i2) else ( i3+1 if i3+1<max(i1,i2) else i3+2 )
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

#####  COMPUTE HOMOGRAPHY MATRIX WITH MINIMUM NUMBER OF CORRESPONDENCES  #####

def compH_tqn(qray,dray,constants):
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
    prm = np.append(np.abs(k)*np.dot(wRd,t),[qYaw,180/np.pi*n])
    if valid: prm = lsqH_tqn(prm,qray,dray,constants)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_tq(qray,dray,constants):
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
    prm = np.append(np.abs(k)*np.dot(wRd,t),qYaw)
    if valid: prm = lsqH_tq(prm,qray,dray,constants)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_tn(qray,dray,constants):
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
    prm = np.append(k*np.dot(wRd,t),180/np.pi*n)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

def compH_t(qray,dray,constants):
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
    prm = np.abs(k)*np.dot(wRd,t)
    valid = valid and geom.vecnorm(prm[:3]) < 5
    return prm, valid

##### SOLVE FOR HOMOGRAPHY MATRIX THAT MINIMIZES REPROJECTION ERROR #####

def lsqH_tqn(prm,qray,dray,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    try: return opt.leastsq(homerrf_tqn,prm,args=(qray,dray,pr,wRd),warning=False)[0]
    except TypeError: return prm

def lsqH_tq(prm,qray,dray,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return opt.leastsq(homerrf_tq,prm,args=(qray,dray,pr,wRd,nYaw),warning=False)[0]

def lsqH_tn(prm,qray,dray,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return opt.leastsq(homerrf_tn,prm,args=(qray,dray,dRq,wRd),warning=False)[0]

def lsqH_t(prm,qray,dray,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return opt.leastsq(homerrf_t,prm,args=(qray,dray,dRq,wRd,nYaw),warning=False)[0]

#####  COMPUTE CORRESPONDENCE ERRORS FOR HOMOGRAPHY COMPUTATIONS  #####

def homerrf_tqn(prm,qray,dray,pr,wRd):
    # set variables
    dRq = np.dot(tp(wRd),geom.RfromYPR(prm[3],pr[0],pr[1]))
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(prm[4]*np.pi/180),0,np.cos(prm[4]*np.pi/180)])
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.append( qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] , qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] )
    return err

def homerrf_tq(prm,qray,dray,pr,wRd,nbear):
    # set variables
    dRq = np.dot(tp(wRd),geom.RfromYPR(prm[3],pr[0],pr[1]))
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(nbear*np.pi/180),0,np.cos(nbear*np.pi/180)])
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.append( qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] , qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] )
    return err

def homerrf_tn(prm,qray,dray,dRq,wRd):
    # set variables
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(prm[3]*np.pi/180),0,np.cos(prm[3]*np.pi/180)])
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.append( qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] , qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] )
    return err

def homerrf_t(prm,qray,dray,dRq,wRd,nbear):
    # set variables
    td = np.dot(tp(wRd),prm[:3])
    nd = -np.dot(tp(wRd),[np.sin(nbear*np.pi/180),0,np.cos(nbear*np.pi/180)])
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(td,nd))
    # Compute homography error
    Hd = tp(np.dot(H,tp(dray)))
    err = np.append( qray[:,0]/qray[:,2]-Hd[:,0]/Hd[:,2] , qray[:,1]/qray[:,2]-Hd[:,1]/Hd[:,2] )
    return err

#####  COMPUTE OVERALL ERROR FOR CORRESPONDENCES  #####

def errH_tqn(prm,qray,dray,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return np.sqrt(np.sum(np.reshape(homerrf_tqn(prm,qray,dray,pr,wRd),[2,-1])**2,0))

def errH_tq(prm,qray,dray,constants):
    Rpr, wRd, qYaw, nYaw = constants
    pr = geom.YPRfromR(Rpr)[1:] # pitch and roll
    return np.sqrt(np.sum(np.reshape(homerrf_tq(prm,qray,dray,pr,wRd,nYaw),[2,-1])**2,0))

def errH_tn(prm,qray,dray,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return np.sqrt(np.sum(np.reshape(homerrf_tn(prm,qray,dray,dRq,wRd),[2,-1])**2,0))

def errH_t(prm,qray,dray,constants):
    wRq, wRd, qYaw, nYaw = constants
    dRq = np.dot(tp(wRd),wRq)
    return np.sqrt(np.sum(np.reshape(homerrf_t(prm,qray,dray,dRq,wRd,nYaw),[2,-1])**2,0))

# computing the normal vector of the planes
def computeNorm(qJpg,dJpg,wRq,wRd,Kq,Kd):
    	# qJpg: query jpeg path
	# dJpg: database jpeg path
	# wRq:  orientation of query in world frame, rotation matrix
	# wRd:  orientation of database in world frame, rotation matrix
	# Kq:   camera matrix for query
	# Kd:   camera matrix for database

    # set parameters
    ymax = .15
    maxerr = .01

    # compute vanishing points
    print 'Extracting query and database vanishing points...'
    qVps, qErrs, qNumi = vanPts.vanishing_points(qJpg,Kq)
    dVps, dErrs, dNumi = vanPts.vanishing_points(dJpg,Kd)

    # move vanishing points to world frame
    qVps = np.dot( wRq , qVps.transpose() ).transpose()
    dVps = np.dot( wRd , dVps.transpose() ).transpose()

    # remove vanishing points with large y components so they're parallel to ground plane
    idx_keep = qVps[:,1] < ymax
    qVps, qErrs, qNumi = qVps[idx_keep,:], qErrs[idx_keep], qNumi[idx_keep]
    idx_keep = dVps[:,1] < ymax
    dVps, dErrs, dNumi = dVps[idx_keep,:], dErrs[idx_keep], dNumi[idx_keep]

    # compute smallest error metric between q/db vanishing point pairs
    # choose the best one
    nbErr = np.inf
    vpQ = np.nan
    vpD = np.nan
    g = np.array([0,1,0])
    for i in xrange(len(qNumi)):
        for j in xrange(len(dNumi)):
            sinQD = np.sqrt( 1 - np.dot(qVps[i,:],dVps[j,:])**2 )
            cosQG = np.abs( np.dot(g,qVps[i,:]) )
            cosDG = np.abs( np.dot(g,dVps[j,:]) )
            err = sinQD * ( cosQG + cosDG )
#            print [sinQD, cosQG, cosDG, err]
#            print err<maxerr
#            print nbErr
#            print err<nbErr
            if err < maxerr and err < nbErr:
                nbErr = err
                vpQ = qVps[i,:]
                vpD = dVps[j,:]

    # point vanishing points the same way and average them
    if np.dot(vpQ,vpD)<0:
        vpD = -vpD
    vp = 1/2. * ( vpQ + vpD )
    vp = vp / alg.norm(vp)

    print vp

    # compute normal vector and face it in the right direction
    n = np.cross( g , [vp[0],0,vp[2]] )
    n = n / alg.norm(n)
    qDir, dDir = np.dot( wRq.transpose() , vp )[2] , np.dot( wRd.transpose() , vp )[2]
    if qDir+dDir < 0:
        n = -n

    # compute bearing of normal vector
    nbear = np.arctan2( n[2] , n[0] )

    return nbear, nbErr

#def essmat_fromR(yd, yq, dR):
#
#    ncand = len(yd)
#    if ncand < 3:
#        raise Exception('Too few correspondences to compute essential matrix.')
#    x, y = np.array(yd), np.array(yq)
#    z = np.array([dot(tp(dR),q) for q in y])
#
#    # RANSAC outer loop
#    infarr = np.array([np.inf,np.inf,np.inf])
#    iter, maxiter = 0, 1000
#    bt, bT, bmask, berr = infarr, infarr, np.zeros(ncand), np.inf
#    nfit, maxerr = 2, 0.05
#    minfit = max(3,int(.05*ncand))
#    while iter < maxiter:
#
#        # fit t with 2 random correspondences
#        idx = np.random.permutation(np.arange(ncand))
#        mx, mz = x[idx[:nfit],:], z[idx[:nfit],:] # nfit random of x,z
#        M = np.zeros([2,3])
#        for i in range(nfit):
#            M[i,:] = np.array( [ mx[i,1]*mz[i,2]-mx[i,2]*mz[i,1] ,
#                            mx[i,2]*mz[i,0]-mx[i,0]*mz[i,2] ,
#                            mx[i,0]*mz[i,1]-mx[i,1]*mz[i,0] ] )
#        u,s,v = alg.svd(M)
#        mt = v[2,:]
#        mT = np.array([[0,-mt[2],mt[1]],
#                  [mt[2],0,-mt[0]],
#                  [-mt[1],mt[0],0]])
#
#        err = np.array([dot(z[i,:],dot(mT,x[i,:])) / (alg.norm(z[i,:])*alg.norm(x[i,:])) for i in range(ncand)])
#        inmask = abs(err) < maxerr
#        ix, iz, im = np.compress(inmask,x,0), np.compress(inmask,z,0), sum(inmask)
#
#        if sum(inmask) >= minfit:
#            M = np.zeros([im,3])
#            for i in range(im):
#                M[i,:] = np.array( [ ix[i,1]*iz[i,2]-ix[i,2]*iz[i,1] ,
#                                ix[i,2]*iz[i,0]-ix[i,0]*iz[i,2] ,
#                                ix[i,0]*iz[i,1]-ix[i,1]*iz[i,0] ] )
#            u,s,v = alg.svd(M)
#            ierr = (s[2]/s[0])**0.5 / im
#            if ierr < berr:
#                bmask, berr = inmask, ierr
#                bt = v[2,:]
#                bT = np.array([[0,-bt[2],bt[1]],
#                          [bt[2],0,-bt[0]],
#                          [-bt[1],bt[0],0]])
#
#        iter += 1
#
#    # essential matrix
#    return (dot(dR,bT), bt, bmask, berr)
#
#
#def TfromR(dR,cf1,cf2):
#    # compute the plane normals
#    npts = len(cf1)
##    print 'Camera frame 1'
##    print cf1
##    print 'Camera frame 1 rotated into cf2'
##    print np.array(dot(cf1,dR))
##    print 'Camera frame 2'
##    print cf2
#    pnorm = [dot(geom.xprodmat(cf1[i,:]),dot(dR,cf2[i,:])) for i in range(npts)]
#    pnorm = np.array([p/alg.norm(p) for p in pnorm]) #normalize
##    print 'Plane normals'
##    print pnorm
#
#    # RANSAC outer loop
#    iter, maxiter = 0, 1000
#    bt, bmask, berr = np.array([np.inf,np.inf,np.inf]), np.ones(npts), np.inf
#    nfit, maxerr = 2, 0.05
#    minfit = max(2,int(.05*npts))
#    while iter < maxiter:
#
#        iter += 1
#
#        # fit t with 2 random correspondences
#        idx = np.random.permutation(np.arange(npts))
#        mpn = pnorm[idx[:nfit],:] # nfit random plane norms
#        if dot(mpn[0,:],mpn[1,:]) > .999:
#            continue
#        mt = dot(geom.xprodmat(mpn[0,:]),mpn[1,:])
#        mt = mt / alg.norm(mt)
#
#        err = dot(pnorm,mt)
#        inmask = abs(err) < maxerr
#        ipn, im = np.compress(inmask,pnorm,0), sum(inmask)
#
#        if sum(inmask) >= minfit:
#            u,s,v = alg.svd(ipn)
#            if len(s) == 2:
#                ierr = np.dot(ipn[0,:],ipn[1,:])
#            else:
#                ierr = s[2] / im
#            if ierr < berr:
#                bmask, berr = inmask, ierr
#                bt = v[2,:]
#
##    # rotation matrix
##    print 'rotation matrix'
##    print dR
##    print 'YPR of rot mat'
##    print 180*np.array(geom.YPRfromR(dR))/np.pi
##
##    # theory translation
##    dt = dot(tp(dR),np.array([0,-2,-10]))
##    dt = dt / alg.norm(dt)
##    print 'theory dt'
##    print dt
##    print 'practice errs'
##    print dot(pnorm,bt)
##    print 'theory errs'
##    print dot(pnorm,dt)
#
#    # translation
#    return (dot(dR,bt), bt, bmask, berr)
#
#
#def scale2view(w3d,q2d,tq,wRq):
#
#    # compute direction of 3d point from query image
#    y = dot(wRq,q2d)
#
#    # drop altitude components
#    x = np.delete(w3d,1)
#    y = np.delete(y,1)
#    z = np.delete(tq,1)
#
#    # set up Aq=k equation representing intersection of 2 lines
#    A = np.array([[y[1],-y[0]],[z[1],-z[0]]])
#    k = y[1]*x[0]-y[0]*x[1]
#    k = np.array([k,0])
#
#    # solve for intersection of 2 lines to get query location
#    q = alg.solve(A,k)
#
#    # compute the corresponding scale factor attached to y
#    s = q[1]/z[1] if q[0]<q[1] else q[0]/z[0]
#
#    return s,q
#
#
##def 3dof_fundmat(Kd, Rd, db2d, Kq, Rq, q2d)
#
#
#
##    return F, R, t
#
#def min_reproj(C, Q, matches, dbimgpath, dbsiftpath):
#
#    # Get image information.
#    info = os.path.join(C.infodir, os.path.basename(dbimgpath)[:-4] + '.info')
#    qsource = render_tags.QueryImageInfo(Q.datasource)
#    dbsource = render_tags.EarthmineImageInfo(dbimgpath, info)
#    map3d = C.pixelmap.open(dbsiftpath)
#
#    # Get 2d pixels and 3d locations of feature inliers
#    matches = [(m['query'][0],m['query'][1],m['db'][0],m['db'][1]) for m in matches]
#    matches = list(set(matches))
#    q2d = [np.array([m[0],m[1]]) for m in matches]
#    db2d = [np.array([m[2],m[3]]) for m in matches]
#    db3d = [map3d[int(d[0]),int(d[1])] for d in db2d]
#    i = 0
#    while i < len(db3d):
#        if db3d[i] is None:
#            q2d.pop(i)
#            db2d.pop(i)
#            db3d.pop(i)
#        else:
#            i = i+1
#    olat,olon,oalt = dbsource.lat,dbsource.lon,dbsource.alt
#    qlat,qlon = qsource.lat,qsource.lon
#    qzx = geom.lltom(olat,olon,qlat,qlon)
#    zx = [geom.lltom(olat,olon,d['lat'],d['lon']) for d in db3d]
#    #print np.array([d['alt'] for d in db3d])
#    #print dbsource.alt
#    y = [dbsource.alt-d['alt'] for d in db3d]
#    xyz = [[zx[i][1],y[i],zx[i][0]] for i in range(len(y))]
#    #print -tp(np.array(y))
#
#    # Set K, e3, Rhat and get list of ci, xi
#    wx,wy = qsource.pgmsize[0], qsource.pgmsize[1]
#    tx,ty = qsource.view_angle[0], qsource.view_angle[1]
#    imsize = np.array([wx,wy])
#    f1, f2 = (wx-1)/(2*np.tan(tx/2)), (wy-1)/(2*np.tan(ty/2))
#    f = (f1+f2) / 2
#    print (wx,wy,tx,f)
#    lsprm = np.concatenate((np.array([f]),imsize),1)
#    K = np.array([[f,0,(wx-1)/2.0],
#             [0,f,(wy-1)/2.0],
#             [0,0,1]])
#    e3 = tp(np.array([[0,0,1]]))
#    y,p,r = qsource.yaw, qsource.pitch, qsource.roll
#    #r = -np.pi/2
#    print [180*y/np.pi,180*p/np.pi,180*r/np.pi]
#    Ry = np.array([[np.cos(y),0,-np.sin(y)],
#              [0,1,0],
#              [np.sin(y),0,np.cos(y)]])
#    Rx = np.array([[1,0,0],
#              [0,np.cos(p),np.sin(p)],
#              [0,-np.sin(p),np.cos(p)]])
#    Rz = np.array([[np.cos(r),np.sin(r),0],
#              [-np.sin(r),np.cos(r),0],
#              [0,0,1]])
#    Rhat = dot(Ry,dot(Rx,Rz))
#    ci = [tp(np.array([[q[0],q[1],1]])) for q in q2d]
#    xi = [tp(np.array([[x[0],x[1],x[2]]])) for x in xyz]
#
#    # Minimize linearized error for T based on orientation R from sensors
#    # Note that T represents translation from camera to global
#    # Global to camera is t = -inv(R)T or T = -Rt
#    # Solution to minimization of linear error is T = -inv(A)*x
#    N = len(xi)
#    A = np.array([[0,0,0],
#             [0,0,0],
#             [0,0,0]])
#    x = tp(np.array([[0,0,0]]))
#    for i in range(N):
#        Mi = K - dot(ci[i],tp(e3))
#        A = A + dot(tp(Mi),Mi)
#        x = x + dot(dot(dot(tp(Mi),Mi),Rhat),xi[i])
#    T = alg.solve(-A,x)
#    t = -dot(tp(Rhat),T)
#    #print [180*p/np.pi,180*y/np.pi,180*r/np.pi]
#    #print np.reshape(t,(3)).tolist()
#
#    # Reshape 3d and 2d points
#    numpts,foo,bar = np.shape(ci)
#    ci = tp(np.reshape(np.array(ci)[:,:2,0],(numpts,-1)))
#    xi = tp(np.reshape(np.array(xi),(numpts,-1)))
#
#    # init orientation
#    rot = np.array([y,p,r])
#
#    # full optimization with least squares
#    errfunc = reproj_errfunc1
#    p0 = np.array([y,p,r,T[0],T[1],T[2]])
#    p0 = np.array([0,0,0])
#    p2, success = opt.leastsq(errfunc,p0,args=(rot,xi,ci,lsprm))
#
#    # optimization of T only with ransac least squares
#    errfunc = reproj_errfunc1
#    p0 = np.array([0,0,0])
#    #rprm=[3,100,10,10,3]
#    minfit,maxerr,stoperr,maxiter,nfit = 50,20,4,100,3
#    p1,rscMask,rscX,rscC,rscErr = ransac_leastsq1(errfunc,p0,rot,xi,ci,lsprm,minfit,maxerr,stoperr,maxiter,nfit)
#    rscX = None
#    tries = 1
#    if rscX is None:
#        minfit,maxerr,stoperr,maxiter,nfit = 20,20,4,50,3
#        p1,rscMask,rscX,rscC,rscErr = ransac_leastsq1(errfunc,p0,rot,xi,ci,lsprm,minfit,maxerr,stoperr,maxiter,nfit)
#        tries = 2
#    if rscX is None:
#        minfit,maxerr,stoperr,maxiter,nfit = 10,20,4,20,3
#        p1,rscMask,rscX,rscC,rscErr = ransac_leastsq1(errfunc,p0,rot,xi,ci,lsprm,minfit,maxerr,stoperr,maxiter,nfit)
#        tries = 3
#    if rscX is None:
#        rscX = xi
#        rscC = ci
#        p1 = p2
#        tries = 0
#
#    # full optimization with ransac least squares
#    errfunc = reproj_errfunc
#    p0 = np.concatenate((rot,p1),1)
#    #rprm=[3,100,10,10,3]
#    minfit,maxerr,stoperr,maxiter,nfit = 50,20,4,100,3
#    p,rscMask,rscX,rscC,rscErr = ransac_leastsq(errfunc,p0,xi,ci,lsprm,minfit,maxerr,stoperr,maxiter,nfit)
#    rscX = None
#    tries = 1
#    if rscX is None:
#        minfit,maxerr,stoperr,maxiter,nfit = 20,20,4,50,3
#        p,rscMask,rscX,rscC,rscErr = ransac_leastsq(errfunc,p0,xi,ci,lsprm,minfit,maxerr,stoperr,maxiter,nfit)
#        tries = 2
#    if rscX is None:
#        minfit,maxerr,stoperr,maxiter,nfit = 10,20,4,20,3
#        p,rscMask,rscX,rscC,rscErr = ransac_leastsq(errfunc,p0,xi,ci,lsprm,minfit,maxerr,stoperr,maxiter,nfit)
#        tries = 3
#    if rscX is None:
#        rscX = xi
#        rscC = ci
#        p = p0
#        tries = 0
#
#    # post-process to get more meaningful return
#    #ry,rx,rz = rot
#    #T = tp(np.array(p[0:3]))
#    ry,rx,rz = p[0:3]
#    T = tp(np.array(p[3:6]))
#    Ry = np.array([[np.cos(ry),0,-np.sin(ry)],
#              [0,1,0],
#              [np.sin(ry),0,np.cos(ry)]])
#    Rx = np.array([[1,0,0],
#              [0,np.cos(rx),np.sin(rx)],
#              [0,-np.sin(rx),np.cos(rx)]])
#    Rz = np.array([[np.cos(rz),np.sin(rz),0],
#              [-np.sin(rz),np.cos(rz),0],
#              [0,0,1]])
#    R = dot(Ry,dot(Rx,Rz))
#    t = -dot(tp(R),T)
#    view = [ry,rx,rz]
#    pos = [t[0],t[1],t[2]]
#    print [180*view[0]/np.pi,180*view[1]/np.pi,180*view[2]/np.pi]
#    #print T.tolist()
#    print pos
#
#    # test position with null space of KR
#    K = np.array([[f,0,(imsize[0]-1)/2.0],
#             [0,f,(imsize[1]-1)/2.0],
#             [0,0,1]])
#    Ry = np.array([[np.cos(ry),0,-np.sin(ry)],
#              [0,1,0],
#              [np.sin(ry),0,np.cos(ry)]])
#    Rx = np.array([[1,0,0],
#              [0,np.cos(rx),np.sin(rx)],
#              [0,-np.sin(rx),np.cos(rx)]])
#    Rz = np.array([[np.cos(rz),np.sin(rz),0],
#              [-np.sin(rz),np.cos(rz),0],
#              [0,0,1]])
#    R = dot(Ry,dot(Rx,Rz))
#    tmp = np.concatenate((R,tp(np.array([T]))),1)
#    MAT = dot(K,np.concatenate((R,tp(np.array([T]))),1))
#    U,S,V = alg.svd(MAT)
#    pt = V[3][:]
#    print [pt[0]/pt[3], p[1]/pt[3], pt[2]/pt[3]]
#
#    # test some 3d to 2d projections
#    param = np.concatenate((p0,np.array([f]),imsize),1)
#    pi = proj3dto2d(param,xi)
#    map = np.array([[ci[0,i],ci[1,i],pi[0,i],pi[1,i]] for i in range(numpts)])
#    err = np.sqrt( (map[:,0]-map[:,2])**2 + (map[:,1]-map[:,3])**2 )
#    avg = np.average(err,0)
#    print 'Average pixel reprojection error for p0: {0}'.format(avg)
#
#    # test some 3d to 2d projections
#    param = np.concatenate((p,np.array([f]),imsize),1)
#    pi = proj3dto2d(param,rscX)
#    map = np.array([[rscC[0,i],rscC[1,i],pi[0,i],pi[1,i]] for i in range(np.shape(rscC)[1])])
#    err = np.sqrt( (map[:,0]-map[:,2])**2 + (map[:,1]-map[:,3])**2 )
#    avg = np.average(err,0)
#    print 'Average pixel reprojection error for inliers: {0}'.format(avg)
#
#    # test some 3d to 2d projections
#    param = np.concatenate((p,np.array([f]),imsize),1)
#    pAll = proj3dto2d(param,xi)
#    mapAll = np.array([[ci[0,i],ci[1,i],pAll[0,i],pAll[1,i]] for i in range(np.shape(ci)[1])])
#    errAll = np.sqrt( (mapAll[:,0]-mapAll[:,2])**2 + (mapAll[:,1]-mapAll[:,3])**2 )
#    avgAll = np.average(errAll,0)
#    print 'Average pixel reprojection error for everything: {0}'.format(avgAll)
#
#    gps_err = ( (t[0]-qzx[1])**2 + (t[2]-qzx[0])**2 )**0.5
#
#    #with open(C.reproj_file,'a') as rp_file:
#    #    print >>rp_file, ''.join([os.path.split(Q.jpgpath)[1][:-4], '\t',
#    #                              str(pos[0]), '\t', str(pos[1]), '\t', str(pos[2]), '\t',
#    #                              str(180*view[0]/np.pi), '\t', str(180*view[1]/np.pi), '\t', str(180*view[2]/np.pi), '\t',
#    #                              str(avg), '\t', str(len(err)), '\t',
#    #                              str(avgAll), '\t', str(len(errAll)), '\t',
#    #                              str(gps_err), '\t', str(tries)])
#
#    return (view,pos,map)
#
#def ransac_leastsq(errf,imodel,x,c,lsprm,minfit=20,maxerr=10,stoperr=2,maxiter=200,nfit=3):
#    # errf : error function for least squares optimization
#    # imodel : initial model parameters R and T
#    # x : 3d locations
#    # c : 2d pixels
#    # lsprm: least squares parameters;
#    #   lsprm[0] = focal length
#    #   lsprm[1] = width of image in pixels
#    #   lsprm[2] = height of image in pixels
#
#    iter = 0
#    bmodel = None
#    bx, bc = None, None
#    bmask = None
#    berr = np.inf
#    nix = 0
#
#    while iter < maxiter and berr > stoperr:
#
#        idx = np.random.permutation(np.arange(np.shape(x)[1]))
#        mx,mc = x[:,idx[:nfit]],c[:,idx[:nfit]] # nfit random of x,c
#        mmodel = opt.leastsq(errf,imodel,args=(mx,mc,lsprm))[0]
#
#        err = errf(mmodel,x,c,lsprm)
#        err = np.reshape(err,[2,len(err)/2])
#        err = np.sqrt( err[0,:]**2 + err[1,:]**2 )
#        inmask = err < maxerr
#        ix,ic = np.compress(inmask,x,1), np.compress(inmask,c,1)
#
#        if np.shape(ix)[1] >= minfit:
#            gmodel = opt.leastsq(errf,mmodel,args=(ix,ic,lsprm))[0]
#            gerr = errf(gmodel,ix,ic,lsprm)
#            gerr = np.reshape(gerr,[2,len(gerr)/2])
#            gerr = np.average(np.sqrt( gerr[0,:]**2 + gerr[1,:]**2 ))
#            if gerr < berr:
#                bmodel = gmodel
#                bx,bc = ix,ic
#                bmask = inmask
#                berr = gerr
#
#        iter += 1
#
#    return (bmodel,bmask,bx,bc,berr)
#
#def reproj_errfunc(p,x,c,lsprm):
#    param = np.concatenate((p,lsprm),1)
#    result = proj3dto2d(param,x) - c
#    result = np.reshape(result,[np.size(result)])
#    return result
#
#def ransac_leastsq1(errf,imodel,rot,x,c,lsprm,minfit=20,maxerr=10,stoperr=2,maxiter=200,nfit=3):
#    # errf : error function for least squares optimization
#    # imodel : initial model parameters R and T
#    # x : 3d locations
#    # c : 2d pixels
#    # lsprm: least squares parameters;
#    #   lsprm[0] = focal length
#    #   lsprm[1] = width of image in pixels
#    #   lsprm[2] = height of image in pixels
#
#    iter = 0
#    bmodel = None
#    bx, bc = None, None
#    bmask = None
#    berr = np.inf
#    nix = 0
#
#    while iter < maxiter and berr > stoperr:
#
#        idx = np.random.permutation(np.arange(np.shape(x)[1]))
#        mx,mc = x[:,idx[:nfit]],c[:,idx[:nfit]] # nfit random of x,c
#        mmodel = opt.leastsq(errf,imodel,args=(rot,mx,mc,lsprm))[0]
#
#        err = errf(mmodel,rot,x,c,lsprm)
#        err = np.reshape(err,[2,len(err)/2])
#        err = np.sqrt( err[0,:]**2 + err[1,:]**2 )
#        inmask = err < maxerr
#        ix,ic = np.compress(inmask,x,1), np.compress(inmask,c,1)
#
#        if np.shape(ix)[1] >= minfit:
#            gmodel = opt.leastsq(errf,mmodel,args=(rot,ix,ic,lsprm))[0]
#            gerr = errf(gmodel,rot,ix,ic,lsprm)
#            gerr = np.reshape(gerr,[2,len(gerr)/2])
#            gerr = np.average(np.sqrt( gerr[0,:]**2 + gerr[1,:]**2 ))
#            if gerr < berr:
#                bmodel = gmodel
#                bx,bc = ix,ic
#                bmask = inmask
#                berr = gerr
#
#        iter += 1
#
#    return (bmodel,bmask,bx,bc,berr)
#
#def reproj_errfunc1(p,rot,x,c,lsprm):
#    param = np.concatenate((rot,p,lsprm),1)
#    result = proj3dto2d(param,x) - c
#    result = np.reshape(result,[np.size(result)])
#    return result
#
#
## Minimize nonlinear error for R,T using least squares regression
#def proj3dto2d(param,xi):
#    ry,rx,rz = param[0:3]
#    T = tp(np.array(param[3:6]))
#    f = param[6]
#    imsize = param[7:]
#    K = np.array([[f,0,(imsize[0]-1)/2.0],
#             [0,f,(imsize[1]-1)/2.0],
#             [0,0,1]])
#    Ry = np.array([[np.cos(ry),0,-np.sin(ry)],
#              [0,1,0],
#              [np.sin(ry),0,np.cos(ry)]])
#    Rx = np.array([[1,0,0],
#              [0,np.cos(rx),np.sin(rx)],
#              [0,-np.sin(rx),np.cos(rx)]])
#    Rz = np.array([[np.cos(rz),np.sin(rz),0],
#              [-np.sin(rz),np.cos(rz),0],
#              [0,0,1]])
#    R = dot(Ry,dot(Rx,Rz))
#    foo,numpts = np.shape(xi)
#    A = np.zeros((3,numpts))
#    p2d = np.zeros((2,numpts))
#    for i in range(numpts):
#        A[:,i] = dot(K,T+dot(R,xi[:,i]))
#        p2d[0,i] = A[0,i]/A[2,i]
#        p2d[1,i] = A[1,i]/A[2,i]
#    return p2d
