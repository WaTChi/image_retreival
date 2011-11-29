# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="aaronh"
__date__ ="$Aug 15, 2011 9:19:51 AM$"

from config import *
import time
import Image
import ImageDraw
import pyflann
import render_tags
import scipy.optimize.minpack as opt
import numpy.random as rand
import reader
import numpy as np
import numpy.linalg as alg
from numpy import array as arr
from numpy import transpose as tp
from numpy import dot
import geom
import os
import pnp
import pickle
import solveHomography
#import vanPts

MAX_PIXEL_DEVIATION = 1
FALLBACK_PIXEL_DEVIATIONS = [2,1]
CONFIDENCE_LEVEL = .99999
BAD_HOMOGRAPHY_DET_THRESHOLD = .005
ROT_THRESHOLD_RADIANS = 0.2 # .1 ~ 5 deg


def highresSift(C, dbmatch, qname):

    # timing info
    start = time.time()

    # set sift paths
    dbsift = os.path.join(C.hiresdir,dbmatch+'sift.txt')
    qsift = os.path.join(C.querydir,'hires',qname+'sift.txt')

    # high res rematch parameters
    maxmatch, maxangle = 10, np.pi/3
    maxratio, maxdist = C.pose_param['maxratio'], C.pose_param['maxdist']

    # rematch file
    filename = qname + ';' + dbmatch + ';maxratio=' + str(maxratio) + \
                    ';maxmatch=' + str(maxmatch) + ';maxdist=' + \
                    str(maxdist/1000) + 'k;maxangle=' + \
                    str(int(round(180/np.pi*maxangle))) + '.npz'
    hrRematchFile = os.path.join( C.querydir, 'hires', 'siftmatch', filename )

    ### HIGH RES REMATCH ###
    matches = {}
    if os.path.isfile(hrRematchFile): # load nearest neighbor data
        print 'Loading high-res sift matches...'
        match_data = np.load(hrRematchFile)
        matches['nmat'] = len(match_data['q2d'])
        matches['numq'] = len(match_data['qidx'])-1
        matches['qidx'] = match_data['qidx']
        matches['q2d'] = match_data['q2d']
        matches['qprm'] = match_data['qprm']
        matches['d2d'] = match_data['d2d']
        matches['dprm'] = match_data['dprm']
        matches['nnd'] = match_data['nnd']
    else: # generate nearest neighbor data and save
        print 'Generating high-res sift matches...'
        rdr = reader.get_reader('sift')
        q = rdr.load_file(qsift)
        db = rdr.load_file(dbsift)
        flann = pyflann.FLANN()
        results, dists = flann.nn(db['vec'], q['vec'], 11, algorithm='linear')
        matches = {'qidx': np.array([]), \
                   'd2d': np.zeros([0,2]), \
                   'dprm': np.zeros([0,2]), \
                   'q2d': np.zeros([0,2]), \
                   'qprm': np.zeros([0,2]), \
                   'nnd': np.array([]) }
        nmat, numq = 0, 0
        for i in xrange(len(results)):
            grads = np.mod(np.array([q[i]['geom'][3]-db[results[i][k]]['geom'][3] for k in range(maxmatch)]),2*np.pi)
            dist = np.float_(np.array(dists[i]))
            high_dist = dist[-1]
            ratios = dist[:maxmatch]/high_dist
            idx = np.nonzero( np.logical_and( ratios<maxratio, \
                np.logical_and( np.array(dists[i][:maxmatch])<maxdist , \
                np.logical_or( grads<maxangle , grads-2*np.pi>-maxangle ) ) ) )[0]
            qadd = len(idx)
            if qadd == 0:
                continue
            q2d = np.tile(np.array([q[i]['geom'][1],q[i]['geom'][0]]),[qadd,1])
            qprm = np.tile(np.array(q[i]['geom'][2:].copy()),[qadd,1])
            d2d = np.array( [ [ db[results[i][k]]['geom'][1] , db[results[i][k]]['geom'][0] ] for k in idx ] )
            dprm = np.array( [ db[results[i][k]]['geom'][2:].copy() for k in idx ] )
            nnd = np.array( [ dists[i][k] for k in idx ] )
            matches['d2d'] = np.append(matches['d2d'],d2d,axis=0)
            matches['dprm'] = np.append(matches['dprm'],dprm,axis=0)
            matches['q2d'] = np.append(matches['q2d'],q2d,axis=0)
            matches['qprm'] = np.append(matches['qprm'],qprm,axis=0)
            matches['nnd'] = np.append(matches['nnd'],nnd,axis=0)
            matches['qidx'] = np.append(matches['qidx'],nmat)
            nmat += qadd
            numq += 1
        matches['qidx'] = np.append(matches['qidx'],nmat)
        matches['nmat'] = nmat
        matches['numq'] = numq
        np.savez( hrRematchFile, qidx=matches['qidx'], d2d=matches['d2d'], dprm=matches['dprm'],
                  q2d=matches['q2d'], qprm=matches['qprm'], nnd=matches['nnd'] )
        # END OF RUNNING NN SEARCH OR LOADING NN FILE
    # fill in other keys
    nmat = matches['nmat']
    numq = matches['numq']
    matches['dray'] = np.zeros([nmat,3])
    matches['qray'] = np.zeros([nmat,3])
    matches['w3d'] = np.zeros([nmat,3])
    matches['plane'] = np.int_( -1 * np.ones(nmat) )
    matches['hmask'] = np.bool_(np.zeros(nmat))
    matches['hvrf'] = 0
    # Print rematch statistics
    print 'Number of query features matched: %d' % numq
    print 'Total number of feature matches: %d' % nmat
    print 'Average number of query matches considered = %.1f' % (float(nmat)/numq)
    print 'High res rematch took %.1f seconds.' % (time.time() - start)
    return matches


def draw_matches(matches, qimg, dbimg, outimg):

    print 'Drawing homography verified sift matches...'

    def scaledown(image, max_height):
        scale = 1.0
        hs = float(image.size[1]) / max_height
        if hs > 1:
            w,h = image.size[0]/hs, image.size[1]/hs
            scale /= hs
            image = image.resize((int(w), int(h)), Image.ANTIALIAS)
        return image, scale

    assert os.path.exists(dbimg)
    assert os.path.exists(qimg)
    a = Image.open(qimg)
    b = Image.open(dbimg)
    if a.mode != 'RGB':
        a = a.convert('RGB')
    if b.mode != 'RGB':
        b = b.convert('RGB')
    height = b.size[1]
    a, scale = scaledown(a, height)
    assert a.mode == 'RGB'
    off = a.size[0]
    target = Image.new('RGBA', (a.size[0] + b.size[0], height))

    def xdrawline((start,stop), color='hsl(20,100%,50%)', off=0):
        start = [start[0] + off, start[1]]
        stop = [stop[0] + off, stop[1]]
        draw.line(start + stop, fill=color, width=1)

    def xdrawcircle((y,x), col='hsl(20,100%,50%)', off=0):
        r = 3
        draw.ellipse((y-r+off, x-r, y+r+off, x+r), outline=col)

    draw = ImageDraw.Draw(target)
    target.paste(a, (0,0))
    target.paste(b, (off,0))

    match_idx = np.nonzero(matches['hmask'])[0] # np.nonzero( np.ones(matches['nmat'] ) )[0] # np.nonzero(matches['hmask'])[0]
    for idx in match_idx:
        start = scale*matches['q2d'][idx,:]
        stop = matches['d2d'][idx,:]
        stop[0] += off
        try:
            xdrawcircle(start,'red')
        except ValueError:
            print start
        xdrawline((start,stop),'orange')
        xdrawcircle(stop,'red')

    target.save(outimg, 'jpeg', quality=90)


def match_info(C, matches, dbmatch, domplane, wRd):

    print 'Filling out query specific information on sift matches...'
    start = time.time()
    nmat = matches['nmat']
    depths = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-depth.png') ) )
    if C.QUERY=='oakland1':
        planes = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-planes.png') ) )
    else:
        planes = -1 * np.ones(depths.shape)
    iremove = np.array([])
    for i in xrange(nmat):
        pix = (int(round(matches['d2d'][i,1])),int(round(matches['d2d'][i,0])))
        depth = depths[pix]
        plane = planes[pix]
        depthfail = depth == 0 or depth == 65535
        planefail = C.pose_param['use_planes'] and plane != domplane
        groundfail = C.pose_param['remove_ground'] and plane == 0
        if depthfail or planefail or groundfail:
            iremove = np.append(iremove,i)
            continue
        matches['plane'][i] = plane if plane != 255 else -1
        depth = depth / 100.
        wray = np.dot( wRd , matches['dray'][i,:] )
        matches['w3d'][i,:] = depth * wray / geom.vecnorm(wray)
    # Remove indices without 3d information
    matches['q2d'] = np.delete(matches['q2d'],iremove,0)
    matches['qprm'] = np.delete(matches['qprm'],iremove,0)
    matches['qray'] = np.delete(matches['qray'],iremove,0)
    matches['d2d'] = np.delete(matches['d2d'],iremove,0)
    matches['dprm'] = np.delete(matches['dprm'],iremove,0)
    matches['dray'] = np.delete(matches['dray'],iremove,0)
    matches['w3d'] = np.delete(matches['w3d'],iremove,0)
    matches['nnd'] = np.delete(matches['nnd'],iremove,0)
    matches['nmat'] -= len(iremove)
    matches['qidx'] = np.sort(np.array(list(set( [ qi-np.sum(iremove<qi) for qi in matches['qidx']] ))))
    matches['numq'] = len(matches['qidx'])-1
    matches['hmask'] = np.bool_(np.zeros(matches['nmat']))
    print 'Removed %d points, leaving %d feature matches' % (len(iremove), matches['nmat'])
    print 'Retrieving 3d points and planes took %.1f seconds.' % (time.time()-start)
    return matches


def getGTpose(C, qname):
    gtFile = open(os.path.join(C.querydir,'gtLatLonNormYaw.txt'),'r')
    gtStr = gtFile.read()
    gtFile.close()
    gtStr = gtStr[gtStr.find(qname):]
    gtStr = gtStr[gtStr.find('\t')+1:]
    idx = gtStr.find('\t')
    qlat = float(gtStr[:idx])
    gtStr = gtStr[idx+1:]
    idx = gtStr.find('\t')
    qlon = float(gtStr[:idx])
    gtStr = gtStr[idx+1:]
    idx = gtStr.find('\t')
    qnorm = float(gtStr[:idx])
    gtStr = gtStr[idx+1:]
    idx = gtStr.find('\n')
    qyaw = float(gtStr[:idx])
    return qlat, qlon, qnorm, qyaw


def find_domplane(C, dbmatch):
    print 'Finding dominant plane...'
    planes = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-planes.png') ) )
    maxplane = np.max(planes[planes!=255])
    plane_count = [ np.sum(np.reshape(planes==k+1,-1)) for k in range(maxplane) ]
    domplane = np.argmax(plane_count) + 1
    second = 0 if maxplane==1 else np.max(np.compress(plane_count!=plane_count[domplane-1],plane_count))
    domsize = plane_count[domplane-1] / np.float(np.size(planes))
    domratio = second / np.float(plane_count[domplane-1])
    return domplane, domsize, domratio


def planefrom3d(C, dbmatch, domplane, Kdinv, wRd):

    print 'Computing dominant plane normal...'

    if domplane == -1: # no dominant plane
        return np.nan, -1

    # get 3d points on plane
    planes = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-planes.png') ) )
    depths = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-depth.png') ) )
    pray = tp( np.dot( wRd , np.dot( Kdinv , tp( np.array( [ [x,y,1] for (y,x) in zip(*np.nonzero(planes==domplane)) ] ) ) ) ) )
    pray = pray / np.tile( tp( [ np.sum( pray**2 , 1 )**0.5 ] ) , [1,3] )
    pdep = tp( np.tile( [ depths[y,x] for (y,x) in zip(*np.nonzero(planes==domplane)) ] , [3,1] ) )
    p3d = np.append( pray * pdep , tp(np.array([np.ones(len(pray))])) )

    # extract least squared normal vector
    n = alg.eig(np.dot(tp(p3d),p3d))[1][2,:3]

    # restrict normal vector to be perpendicular to gravity
    g = np.array([0,1,0])
    n = geom.normalrows( n - np.inner(n,g) * g )

    # get plane distance and flip n if needed
    pd = np.mean(np.dot(p3d[:,:3],n))
    if pd<0:
        n, pd = -n, -pd

    return n, pd


def testplanefrom3d(p3d):

    print 'Computing dominant plane normal...'

    # extract least squared normal vector
    u,s,vt = alg.svd(p3d)
    n = geom.normalrows(vt[2,:3])

    # restrict normal vector to be perpendicular to gravity
    g = np.array([0,1,0])
    n = geom.normalrows( n - np.inner(n,g) * g )

    # get plane distance and flip n if needed
    pd = np.mean(np.dot(p3d[:,:3],n))
    if pd<0:
        n, pd = -n, -pd

    return n, pd


def estimate_pose(C, Q, dbmatch, gtStatus=None):

    # get pose parameters
    param = C.pose_param

    # get high res db image and sift paths
    dbinfo = os.path.join(C.hiresdir, dbmatch + '.info')
    dbimg = os.path.join(C.hiresdir,dbmatch+'.jpg')
    dbsift = os.path.join(C.hiresdir,dbmatch+'sift.txt')
    dbsource = render_tags.EarthmineImageInfo(dbimg, dbinfo)

    # get high res query information
    qname = Q.name
    qimg = os.path.join(C.querydir,'hires',qname+'.jpg')
    qsift = os.path.join(C.querydir,'hires',qname+'sift.txt')
    qsource = render_tags.QueryImageInfo(Q.datasource)

    # get high res sift rematch
    matches = highresSift(C, dbmatch, qname)

    # temporarily get qyaw to see performance with "good" yaw
    qlat, qlon, qnorm, qyaw = getGTpose(C, qname)

    # Set Kq, Rq
    wx,wy = qsource.image.size
    fov = qsource.view_angle[0]
    Kq = geom.cameramat(wx, wy, fov)
    Kqinv = alg.inv(Kq)
    y,p,r = qsource.yaw, qsource.pitch, qsource.roll
    yaw_used = qyaw
    Rq = geom.RfromYPR(yaw_used,p,r) # camera orientation (camera to world)

    # temporarily get yaw error
    yaw_error = int(abs(np.mod(qyaw-yaw_used,360)))
    yaw_error = yaw_error if yaw_error<180 else abs(yaw_error-360)

    # Set Kd, Rd, and db position
    wx,wy = dbsource.image.size
    fov = dbsource.fov
    Kd = geom.cameramat(wx, wy, fov)
    Kdinv = alg.inv(Kd)
    y,p,r = dbsource.yaw, dbsource.pitch, dbsource.roll
    Rd = geom.RfromYPR(y,p,r) # db camera orientation (camera to world)
    olat,olon,oalt = dbsource.lat,dbsource.lon,dbsource.alt # database location

    # Check for dominant plane and set flags
    method = param['mode']
    if method=='ess' and param['use_planes']:
        print 'Cannot use plane information with Essential matrix.'
        param['use_planes'] = False
    if method == 'comb' and not param['use_planes']:
        print 'Muse use plane information with a combination approach.'
        param['use_planes'] = True
    domplane, domsize, domratio = find_domplane(C, dbmatch) if param['use_planes'] else (-1, 0, 1)
    if domsize < param['min_domsize']:
        print 'Dominant plane not large enough to restrict features.'
        param['use_planes'] = False
    if method=='comb':
        method = 'hom' if param['use_planes'] else 'ess'

    # Fill out match information
    nmat = matches['nmat']
    matches['qray'] = tp(dot(Kqinv,np.append(tp(matches['q2d']),[np.ones(nmat)],0)))
    matches['dray'] = tp(dot(Kdinv,np.append(tp(matches['d2d']),[np.ones(nmat)],0)))
    match_info(C, matches, dbmatch, domplane, Rd)

    # approximate normal vector from 3d points and compute error
    norm_err = 180
#    if C.QUERY == 'oakland1':
#        print 'Dominant plane: %d' % domplane
#        n3d, pd3d = planefrom3d(C, dbmatch, domplane, Kdinv, Rd)
#        nbear3d = 180 / np.pi * np.arctan2(n3d[0],n3d[2])
#        norm_err = int(abs(np.mod(180+qnorm-nbear3d,360)))
#        norm_err = norm_err if norm_err<180 else abs(norm_err-360)
#        print 'Plane normal error from 3d points: %d degrees' % norm_err

    # Get estimated ground truth query location and normal direction
    qlat, qlon, qnorm, qyaw = getGTpose(C, qname)
    qzx = geom.lltom(olat,olon,qlat,qlon)

    # Retrieve GPS query location
    glat,glon = qsource.lat,qsource.lon
    gzx = geom.lltom(olat,olon,glat,glon)

#    # Solve for normal vector using vanishing points
#    if vpN:
#        vbear, vbErr = computeNorm(Q.jpgpath,hrdbimg,Rq,Rd,Kq,Kd)
#    else:
#        vbear = np.nan
#    vbeardeg = 180/np.pi*vbear - 180*round((vbear*180/np.pi-qbear)/180.)
#    vAng_err = np.nan if np.isnan(vbear) else int(abs(vbeardeg-qbear))
#    print vbear
#    print vAng_err

    # Solve for query pose using constrained homography
    matches, ct, nbear = constrainedHomography(matches,Rd,Rq,nbear=np.nan,maxerr=param['inlier_error'],maxiter=param['ransac_iter'])
    if matches['hvrf']==0 or alg.norm(ct)>50:
        print 'Attempting to solve constrained homography with relaxed parameters.'
        matches, ct, nbear = constrainedHomography(matches,Rd,Rq,nbear=np.nan,maxerr=2*param['inlier_error'],maxiter=param['ransac_iter'])
    if matches['hvrf']==0 or alg.norm(ct)>50:
        print 'Attempting to solve constrained homography with even more relaxed parameters.'
        matches, ct, nbear = constrainedHomography(matches,Rd,Rq,nbear=np.nan,maxerr=5*param['inlier_error'],maxiter=param['ransac_iter'])
    nbear = np.mod( 180 + 180/np.pi * nbear , 360 )

    # compute location errors wrt estimated query locations
    ch_err = ( (ct[0]-qzx[1])**2 + (ct[2]-qzx[0])**2 )**0.5
    gps_err = ( (gzx[1]-qzx[1])**2 + (gzx[0]-qzx[0])**2 )**0.5

    # compute the angle difference between T and ground truth translation
    tAng_err = int(abs( 180/np.pi * np.arccos( (ct[0]*qzx[1]+ct[2]*qzx[0]) / (alg.norm([ct[0],ct[2]])*alg.norm(qzx)) ) ))

    # compute the plane normal angle error
    nAng_err = np.nan if np.isnan(nbear) or np.isnan(qnorm) else int(abs(nbear-qnorm))
#    vAng_err = np.nan if np.isnan(vbear) else int(abs(vbear-qbear))

    # write pose estimation results to file
    with open(param['pose_file'],'a') as pose_file:
        print >>pose_file, '\t'.join([qname, str(ch_err), str(gps_err), \
            str(tAng_err), str(nAng_err), str(yaw_error), str(matches['hvrf']), \
            str(matches['numq']), str(matches['nmat']), str(matches['herr']), \
            str(ct[0]), str(ct[2]), str(norm_err)])
#        print >>pose_file, '\t'.join([qname, str(ch_err), str(gps_err), \
#            str(tAng_err), str(nAng_err), str(matches['hvrf']), \
#            str(matches['numq']), str(matches['nmat']), str(matches['herr']), \
#            str(ct[0]), str(ct[2]), str(geom.YPRfromR(Rq))])

    # draw matches
    close = ch_err < 10
    imgpath = os.path.join( param['resultsdir'] , qname + ';gt' + str(gtStatus) + ';' + 'close' + str(close) + ';err=' + str(ch_err) + ';tAng=' + str(tAng_err) + ';nAng=' + str(nAng_err) + ';yAng=' + str(yaw_error) + ';' + dbmatch + '.jpg')
    draw_matches(matches, qimg, dbimg, imgpath)

    print 'YPR of wRq: ' + str(geom.YPRfromR(Rq))
    print 'YPR of wRd: ' + str(geom.YPRfromR(Rd))
    print 'YPR of dRq: ' + str(geom.YPRfromR(np.dot(tp(Rd),Rq)))

    input = (Rq, Rd)


    return ct, ch_err, matches, input

def constrainedHomography(matches, wRd, wRq, nbear=np.nan, maxerr=.05, maxiter=1000):

    print 'Solving constrained homography...'
    start = time.time()

    # Set variables
    nmat, numq = matches['nmat'], matches['numq']
    dRq = dot(tp(wRd),wRq)
    infarr3 = arr([np.inf,np.inf,np.inf])
    T, bp, bmask, bnumi, berr, nbear = infarr3, np.zeros(5), np.zeros(nmat), 0, np.inf, np.nan

    # Set Ransac parameters
    minsep = 0.1 # minimum separation between 2 random RANSAC points
    maxminfit = 40 # maximum threshold for minimum fit
    minfit = min( maxminfit , max( 3 , int(matches['numq']**0.5) ) )
    iter, stoperr = 0, .01*maxerr
    spd = 0 # 1e-2 # scale for plane distance errors relative to homography reprojection errors
    knownN = not np.isnan(nbear)

    # Ransac loop to eliminate outliers with homography
    # Solves homography matrix for homography matrix H=qRd(I+rn') using y ~ Hx
    qray, dray, w3d, qidx = matches['qray'], matches['dray'], matches['w3d'], matches['qidx']
    while iter < maxiter and bnumi < np.inf:
        iter += 1
        # nfit random correspondences
        i0 = rand.randint(0,nmat)
        i1 = rand.randint(0,nmat-1)
        i1 = i1+1 if i1>=i0 else i1
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
            if geom.vecnorm(mp[:3]) > 3:
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
#            mp = opt.leastsq(homerrf_RNP,mp,args=([mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),warning=False)[0]
#            verr = np.reshape(homerrf_RNP(mp,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
#            if (abs(verr[:,:2]>stoperr)).any() or (abs(verr[:,2])>spd).any(): # verify homography solution and 3d world points
#                continue
            if (abs(verr[0,2]-verr[1,2])>spd).any():
                continue
            if geom.vecnorm(mp[:3]) > 3:
                continue
            errs = np.sum(np.reshape(homerrf_RNP(mp,matches['qray'],matches['dray'],matches['w3d'],dRq,wRd,spd),[-1,3])**2,1)**0.5
        emask = errs < maxerr
        imask = np.bool_(np.zeros(matches['nmat']))
        for i in xrange(numq):
            dmask = emask[qidx[i]:qidx[i+1]]
            if np.sum(dmask) == 0:
                continue
            derrs = errs[qidx[i]:qidx[i+1]]
            imask[qidx[i]+np.argmin(derrs)] = True
#        imask = emask
        numi = np.sum(imask)
        if numi > bnumi:
            bp = mp
            bmask = imask
            bnumi = numi
    # Guided matching
    numi, imask, mp = bnumi, bmask, bp
    last_numi, gm_iter = 0, 0
    while last_numi != numi and numi >= minfit:
        if gm_iter >= 10 and last_numi >= numi:
            break
        last_numi = numi
        gm_iter += 1
        iq = qray[imask,:]
        id = dray[imask,:]
        iw = w3d[imask,:]
        try:
            mp = opt.leastsq(homerrf_RNP,mp,args=(iq,id,iw,dRq,wRd,spd),warning=False)[0]
        except TypeError:
            print numi
            print last_numi
            print gm_iter
            print iq
            print mp
            raise TypeError
        errs = np.sum(np.reshape(homerrf_RNP(mp,matches['qray'],matches['dray'],matches['w3d'],dRq,wRd,spd),[-1,3])**2,1)**0.5
        emask = errs < maxerr
        imask = np.bool_(np.zeros(matches['nmat']))
        for i in xrange(numq):
            dmask = emask[qidx[i]:qidx[i+1]]
            if np.sum(dmask) == 0:
                continue
            derrs = errs[qidx[i]:qidx[i+1]]
            imask[qidx[i]+np.argmin(derrs)] = True
        numi = np.sum(imask)
#    if numi >= minfit and alg.norm(homerrf_RNP(mp,iq,id,iw,dRq,wRd,spd)) * numq / numi < berr:
    iq = qray[imask,:]
    id = dray[imask,:]
    iw = w3d[imask,:]
    berr = alg.norm(homerrf_RNP(mp,iq,id,iw,dRq,wRd,spd)) * numq / numi
    bmask = imask
    bnumi = numi
    bp = mp
    T = mp[:3] # mp[4]*mp[:3]
    nbear = bp[3]
    matches['Hprm'] = bp
    matches['Hmat'] = dot(tp(dRq),np.eye(3,3)-dot(r,tp(n)))
    matches['Tvec'] = T
    matches['Pnorm'] = arr([np.sin(bp[3]),0,np.cos(bp[3])])
    matches['hmask'] = bmask
    matches['herr'] = berr
    matches['hvrf'] = sum(bmask)
    if matches['hvrf'] == 0:
        print 'Constrained homography failed.'
    else:
        print 'Result from error metric choosing best inlier set: %f' % matches['herr']
        print 'Number of inliers / total correspondences: ' + str(matches['hvrf']) + ' / ' + str(matches['nmat'])
    print 'Constrained homography took %.1f seconds.' % (time.time()-start)

    return matches, T, nbear


def homerrf_RP(prm,nbear,q,d,w,dRq,wRd,spd):

    # Construct homography matrix from parameters
    r = dot(tp(wRd),prm[:3])
    wn = arr([np.sin(nbear),0,np.cos(nbear)])
    n = dot(tp(wRd),wn)
    H = dot(tp(dRq),np.eye(3,3)-np.outer(r,n))
    # Compute homography error
    q = np.reshape(arr(q),[-1,3])
    d = np.reshape(arr(d),[-1,3])
    w = np.reshape(arr(w),[-1,3])
    Hd = tp(dot(H,tp(d)))
    err = arr( [ [ q[i,0]/q[i,2]-Hd[i,0]/Hd[i,2] , q[i,1]/q[i,2]-Hd[i,1]/Hd[i,2] , spd*(prm[3]-dot(wn,w[0,:])) ] for i in range(len(q)) ] )
    return np.reshape(err,3*len(q))


def homerrf_RNP(prm,q,d,w,dRq,wRd,spd):

    # Construct homography matrix from parameters
    r = dot(tp(wRd),prm[:3])
    wn = arr([np.sin(prm[3]),0,np.cos(prm[3])])
    n = dot(tp(wRd),wn)
    H = dot(tp(dRq),np.eye(3,3)-np.outer(r,n))
    # Compute homography error
    q = np.reshape(arr(q),[-1,3])
    d = np.reshape(arr(d),[-1,3])
    w = np.reshape(arr(w),[-1,3])
    Hd = tp(dot(H,tp(d)))
    err = arr( [ [ q[i,0]/q[i,2]-Hd[i,0]/Hd[i,2] , q[i,1]/q[i,2]-Hd[i,1]/Hd[i,2] , spd*(prm[4]-dot(wn,w[0,:])) ] for i in range(len(q)) ] )
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
#    zerof = lambda e: (1/m[0])*np.cos(e-f[0]) - (1/m[1])*np.cos(e-f[1])
#    e = opt.fsolve(zerof,(f[0]+f[1])/2)
#    e = np.mod(e,2*np.pi)
#    k = m[0] / np.cos(e-f[0])
    errf = lambda prm,argm,argf: prm[0]-argm/np.cos(prm[1]-argf)
    ke_init = arr([1.2*np.mean(m),np.mean(f)])
    k, e = tuple( opt.leastsq(errf,ke_init,args=(m,f))[0] )
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
