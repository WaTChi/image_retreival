#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="Aaron"
__date__ ="$Nov 3, 2011 9:04:33 AM$"

import numpy as np
import scipy as sp
import scipy.optimize.minpack as opt
import numpy.linalg as alg
from numpy import transpose as tp
import geom

#from config import *
#import time
#import Image
#import ImageDraw
#import pyflann
#import render_tags
#import scipy.optimize.minpack as opt
#import numpy.random as rand
#import reader
#import numpy as np
#import numpy.linalg as alg
#from numpy import array as arr
#from numpy import transpose as tp
#from numpy import dot
#import geom
#import os
#import pnp
#import pickle
##import vanPts


def constrainedHomography(matches, wRd, wRq, nbear=np.nan, maxerr=.05, maxiter=1000):

    print 'Solving constrained homography...'
    start = time.time()

    # Set variables
    rngN = range(matches['nmat'])
    dRq = dot(tp(wRd),wRq)
    infarr3 = arr([np.inf,np.inf,np.inf])
    T, berr, nbear = infarr3, np.inf, np.nan

    # Set Ransac parameters
    minsep = 0.1 # minimum separation between 2 random RANSAC points
    maxminfit = 40 # maximum threshold for minimum fit
    minfit = min( maxminfit , max( 3 , int(matches['numq']**0.5) ) )
    iter, stoperr = 0, .01*maxerr
    spd = 1e-3 # scale for plane distance errors relative to homography reprojection errors
    knownN = not np.isnan(nbear)
    print nbear
    print knownN

    # Ransac loop to eliminate outliers with homography
    # Solves homography matrix for homography matrix H=qRd(I+rn') using y ~ Hx
    qray, dray, w3d, qidx = matches['qray'], matches['dray'], matches['w3d'], matches['qidx']
    while iter < maxiter:
        iter += 1
        # nfit random correspondences
        i0, i1 = tuple(rand.permutation(rngN)[:2])
        mq1, mq2 = matches['qray'][i0,:], matches['qray'][i1,:]
        md1, md2 = matches['dray'][i0,:], matches['dray'][i1,:]
        mw1, mw2 = matches['w3d'][i0,:], matches['w3d'][i1,:]
        mq_dist = alg.norm(mq1-mq2)
        if mq_dist < minsep:
            continue
        if knownN:
            r, k, n = lsqHprm_R([mq1,mq2],[md1,md2],dRq,wRd,nbear) # directly solve to get initial conditions
            pd = np.mean([dot(n,mw1),dot(n,mw2)])
            mp = np.append(dot(wRd,k*r),pd) # initial conditions for reprojection error minimization
            verr = np.reshape(homerrf_RP(mp,nbear,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2])>stoperr).any():
                mp[:3] = dot(wRd,-k*r)
                verr = np.reshape(homerrf_RP(mp,nbear,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2])>stoperr).any():
                continue
            mp = opt.leastsq(homerrf,mp,args=([mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),warning=False)[0]
            verr = np.reshape(homerrf(mp,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2]>stoperr)).any() or (abs(verr[:,2])>spd).any(): # verify homography solution and 3d world points
                continue
            errs = np.sum(np.reshape(hom_errf(mp,matches['qray'],matches['dray'],matches['w3d'],dRq,wRd,spd),[-1,3])**2,1)**0.5
        else:
            r, k, e, n = lsqHprm_RN([mq1,mq2],[md1,md2],dRq,wRd) # directly solve to get initial conditions
            pd = np.mean([dot(n,mw1),dot(n,mw2)])
            mp = np.append(dot(wRd,k*r),[e,pd]) # initial conditions for reprojection error minimization
            verr = np.reshape(homerrf_RNP(mp,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2])>stoperr).any():
                mp[:3] = dot(wRd,-k*r)
                verr = np.reshape(homerrf_RNP(mp,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2])>stoperr).any():
                continue
            mp = opt.leastsq(homerrf_RNP,mp,args=([mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),warning=False)[0]
            verr = np.reshape(homerrf_RNP(mp,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2]>stoperr)).any() or (abs(verr[:,2])>spd).any(): # verify homography solution and 3d world points
                continue
            errs = np.sum(np.reshape(homerrf_RNP(mp,matches['qray'],matches['dray'],matches['w3d'],dRq,wRd,spd),[-1,3])**2,1)**0.5
        emask = errs < maxerr
        imask = np.bool_(np.zeros(matches['nmat']))
        print len(emask)
        print qidx
        for i in rngN:
            qmask = emask[qidx[i]:qidx[i+1]]
            if np.sum(qmask) == 0:
                continue
            qerrs = errs[qidx[i]:qidx[i+1]]
            imask[qidx[i]+np.argmin(qerrs)] = True
        numi = np.sum(imask)
        if numi >= minfit:
            iq = qray[imask,:]
            id = dray[imask,:]
            iw = w3d[imask,:]
            ip = opt.leastsq(homerrf_RNP,mp,args=(iq,id,iw,dRq,wRd,spd),warning=False)[0]
            ierr = alg.norm(homerrf_RNP(ip,iq,id,iw,dRq,wRd,spd)) / numi
            if ierr < berr:
                berr = ierr
                T = ip[4]*ip[:3]
                nbear = ip[3]
                matches['Hprm'] = ip
                matches['Hmat'] = dot(tp(dRq),np.eye(3,3)-dot(r,tp(n)))
                matches['Tvec'] = T
                matches['Pnorm'] = arr([np.sin(ip[3]),0,np.cos(ip[3])])
                matches['hmask'] = imask
                matches['herr'][imask,:] = np.sum(np.reshape(homerrf_RNP(ip,iq,id,iw,dRq,wRd,spd),[-1,3])**2,1)**0.5
                matches['hvrf'] = sum(imask)
        if berr < stoperr:
            break
    if matches['hvrf'] == 0:
        print 'Constrained homography failed.'
    else:
        print 'Result from error metric choosing best inlier set: %f' % berr
        print 'Average reprojection and 3d plane error: ' + str(np.mean(matches['herr'][matches['hmask']]))
        print 'Number of inliers / total correspondences: ' + str(matches['hvrf']) + ' / ' + str(matches['nmat'])
    print 'Constrained homography took %.1f seconds.' % (time.time()-start)

    return matches, T, nbear

def lsqHom(qray,dray,w3d,dRq,wRd):
    pass

def hom2pt(q,d,dRq,wRd,nbear=np.nan):

    # Set variables
    rng = [0,1]
    x = [d[i]/alg.norm(d[i]) for i in rng]
    y = [q[i]/alg.norm(q[i]) for i in rng]
    a = [dot(dRq,y[i]) for i in rng]
    b = [dot(wRd,x[i]) for i in rng]
    c = [np.cross(a[i],x[i]) for i in rng]

    # Compute direction of translation
    r = np.cross(c[0],c[1])
    r = r / alg.norm(r) # direction of translation

    # Compute magnitude of translation and bearing of normal vector
    if np.isnan(nbear): # normal vector unknown:
        m = [alg.norm(c[i])/(alg.norm(np.cross(a[i],r))*alg.norm(b[i][[0,2]])) for i in rng]
        f = [np.arctan2(b[i][0],b[i][2]) for i in rng]
        zerof = lambda e: (1/m[0])*np.cos(e-f[0]) - (1/m[1])*np.cos(e-f[1])
        e = opt.fsolve(zerof,(f[0]+f[1])/2)
        e = np.mod(e,2*np.pi)
        k = m[0] / np.cos(e-f[0])
        fe = np.mod(e-np.mean(f),2*np.pi)
    else: # normal vector known
        pass

    if fe > np.pi/2 and fe < 3*np.pi/2:
        e = np.mod(e+np.pi,2*np.pi)
    if np.mean(arr([dot(r,a[i]-x[i]) for i in rng_npts])) < 0:
        r = -r
    n = arr([np.sin(e), 0, np.cos(e)])

    return r, k, e, n

def homerrf(param, qray, dray, w3d, dRq, wRd, scale_pd=0, nbear=np.nan):

    # Construct homography matrix from parameters
    hom_translate = np.dot(tp(wRd),param[:3]) # translation in db frame
    nbear, plane_distance = param[3], param[4] if len(param)==5 else nbear, param[3]
    world_normal = np.array([np.sin(nbear),0,np.cos(nbear)]) # world normal
    normal = np.dot(tp(wRd),world_normal) # normal in db frame
    H = np.dot(tp(dRq),np.eye(3,3)-np.outer(hom_translate,normal)) # homography matrix

    # Compute homography error
    Hd = tp(dot(H,tp(dray))) # mapped database points using homography matrix
    err = np.zeros([len(q),3])
    err[:,0] = qray[:,0]-Hd[:,0]/Hd[:,2]
    err[:,1] = qray[:,1]-Hd[:,1]/Hd[:,2]
    err[:,2] = scale_pd * ( plane_distance-np.dot(w3d,world_normal) )
    return np.reshape(err,3*len(q))


def homerrf_RP(prm,nbear,q,d,w,dRq,wRd,spd):

    # Construct homography matrix from parameters
    r = np.dot(tp(wRd),prm[:3])
    wn = arr([np.sin(nbear),0,np.cos(nbear)])
    n = dot(tp(wRd),wn)
    H = dot(tp(dRq),np.eye(3,3)-np.outer(r,n))
    # Compute homography error
    q = np.reshape(arr(q),[-1,3])
    d = np.reshape(arr(d),[-1,3])
    w = np.reshape(arr(w),[-1,3])
    Hd = tp(dot(H,tp(d)))
    err = arr( [ [ q[i,0]-Hd[i,0]/Hd[i,2] , q[i,1]-Hd[i,1]/Hd[i,2] , spd*(prm[3]-dot(wn,w[0,:])) ] for i in range(len(q)) ] )
    return np.reshape(err,3*len(q))


def homerrf_RNP(prm,q,d,w,dRq,wRd,spd):

    # Construct homography matrix from parameters
    r = np.dot(tp(wRd),prm[:3])
    wn = arr([np.sin(prm[3]),0,np.cos(prm[3])])
    n = dot(tp(wRd),wn)
    H = dot(tp(dRq),np.eye(3,3)-np.outer(r,n))
    # Compute homography error
    q = np.reshape(arr(q),[-1,3])
    d = np.reshape(arr(d),[-1,3])
    w = np.reshape(arr(w),[-1,3])
    Hd = tp(dot(H,tp(d)))
    err = arr( [ [ q[i,0]-Hd[i,0]/Hd[i,2] , q[i,1]-Hd[i,1]/Hd[i,2] , spd*(prm[4]-dot(wn,w[0,:])) ] for i in range(len(q)) ] )
    return np.reshape(err,3*len(q))


def lsqHprm_RN(q,d,dRq,wRd):

    # Set variables
    npts = len(q)
    rng_npts = range(npts)
    x = [d[i]/alg.norm(d[i]) for i in rng_npts]
    y = [q[i]/alg.norm(q[i]) for i in rng_npts]
    a = [dot(dRq,y[i]) for i in rng_npts]
    b = [dot(wRd,x[i]) for i in rng_npts]
    c = [np.cross(a[i],x[i]) for i in rng_npts]

    # Compute homography parameters
    r = np.cross(c[0],c[1])
    r = r / alg.norm(r) # direction of translation
    m = [alg.norm(c[i])/(alg.norm(np.cross(a[i],r))*alg.norm(b[i][[0,2]])) for i in rng_npts]
    f = [np.arctan2(b[i][0],b[i][2]) for i in rng_npts]
    zerof = lambda e: (1/m[0])*np.cos(e-f[0]) - (1/m[1])*np.cos(e-f[1])
    e = opt.fsolve(zerof,(f[0]+f[1])/2)
    e = np.mod(e,2*np.pi)
    k = m[0] / np.cos(e-f[0])
#    errf = lambda prm,argm,argf: prm[0]-argm/np.cos(prm[1]-argf)
#    ke_init = arr([1.2*np.mean(m),np.mean(f)])
#    k, e = tuple( opt.leastsq(errf,ke_init,args=(m,f))[0] )
    fe = np.mod(e-np.mean(f),2*np.pi)
    if fe > np.pi/2 and fe < 3*np.pi/2:
        e = np.mod(e+np.pi,2*np.pi)
    if np.mean(arr([dot(r,a[i]-x[i]) for i in rng_npts])) < 0:
        r = -r
    n = arr([np.sin(e), 0, np.cos(e)])

    return r, k, e, n


def lsqHprm_R(q,d,dRq,wRd,nbear):

    # Set variables
    npts = len(q)
    rng_npts = range(npts)
    x = [d[i]/alg.norm(d[i]) for i in rng_npts]
    y = [q[i]/alg.norm(q[i]) for i in rng_npts]
    a = [dot(dRq,y[i]) for i in rng_npts]
    b = [dot(wRd,x[i]) for i in rng_npts]
    c = [np.cross(a[i],x[i]) for i in rng_npts]

    # Compute homography parameters
    r = np.cross(c[0],c[1])
    r = r / alg.norm(r) # direction of translation
    m = np.array([alg.norm(c[i])/(alg.norm(np.cross(a[i],r))*alg.norm(b[i][[0,2]])) for i in rng_npts])
    f = np.array([np.arctan2(b[i][0],b[i][2]) for i in rng_npts])
    k = np.mean( m / np.cos(nbear-f) )
    if np.mean(arr([dot(r,a[i]-x[i]) for i in rng_npts])) < 0:
        r = -r
    n = arr([np.sin(nbear), 0, np.cos(nbear)])

    return r, k, n

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