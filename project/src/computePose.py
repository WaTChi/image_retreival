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
import reader
import numpy as np
import numpy.linalg as alg
from numpy import transpose as tp
import geom
import os
import solveHomography
import solveEssMatrix
import vp_analysis

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
    matches['imask'] = np.bool_(np.zeros(nmat))
    matches['hvrf'] = 0
    # Print rematch statistics
    print 'Number of query features matched: %d' % numq
    print 'Total number of feature matches: %d' % nmat
    print 'Average number of database matches considered = %.1f' % (float(nmat)/numq)
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

    match_idx = np.nonzero(matches['imask'])[0] # np.nonzero( np.ones(matches['nmat'] ) )[0]
    for idx in match_idx:
        start = scale*matches['q2d'][idx,:]
        stop = matches['d2d'][idx,:]
        stop[0] += off
        xdrawcircle(start,'red')
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
        matches['w3d'][i,:] = depth * wray
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
    matches['imask'] = np.bool_(np.zeros(matches['nmat']))
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
    h,w = planes.shape
    centerplane = planes[int(h/2),int(w/2)]
    maxplane = np.max(planes[planes!=255])
    plane_count = [ np.sum(np.reshape(planes==k+1,-1)) for k in range(maxplane) ]
    domplane = np.argmax(plane_count) + 1
    if domplane != centerplane:
        return -1, 0, 1
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
    p3d = np.append( pray * pdep , tp(np.array([np.ones(len(pray))])) , 1 )

    # extract least squared normal vector
    try:
        neig = alg.eig(np.dot(tp(p3d),p3d))
        n = neig[1][:3,np.argmin(neig[0])]
    except np.linalg.linalg.LinAlgError:
        print p3d.shape
        print np.dot(tp(p3d),p3d)
        neig = alg.eig(np.dot(tp(p3d),p3d))
        n = neig[1][:3,np.argmin(neig[0])]
    #n = alg.eig(np.dot(tp(p3d),p3d))[1][2,:3]

    # restrict normal vector to be perpendicular to gravity
    g = np.array([0,1,0])
    n = geom.normalrows( n - np.inner(n,g) * g )

    # get plane distance and flip n if needed
    pd = np.mean(np.inner(p3d[:,:3],n))
    if pd>0:
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

    # temporarily get tyaw to see performance with "good" yaw
    tlat, tlon, tnorm, tyaw = getGTpose(C, qname) # true lat, lon, normal, yaw

    # Set Kq, Rq
    wx,wy = qsource.image.size
    fov = qsource.fov
    Kq = geom.cameramat(wx, wy, fov)
    Kqinv = alg.inv(Kq)
    cyaw,p,r = qsource.yaw, qsource.pitch, qsource.roll # cyaw - cell phone yaw
    #cyaw=tyaw
    Rq = geom.RfromYPR(cyaw,p,r) # camera orientation (camera to world)

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
    matches['qray'] = geom.normalrows(tp(np.dot(Kqinv,np.append(tp(matches['q2d']),[np.ones(nmat)],0))))
    matches['dray'] = geom.normalrows(tp(np.dot(Kdinv,np.append(tp(matches['d2d']),[np.ones(nmat)],0))))
    matches = match_info(C, matches, dbmatch, domplane, Rd)

    # approximate normal vector from 3d points and compute error
    norm_err = -1
    print 'Dominant plane: %d' % domplane
    if C.QUERY == 'oakland1' and param['use_planes']:
        n3d, pd3d = planefrom3d(C, dbmatch, domplane, Kdinv, Rd)
        nbear3d = 180 / np.pi * np.arctan2(n3d[0],n3d[2])
        norm_err = int(abs(np.mod(tnorm-nbear3d,360)))
        norm_err = norm_err if norm_err<180 else abs(norm_err-360)
        print 'Plane normal error from 3d points: %d degrees' % norm_err

    # Get estimated ground truth query location and normal direction
    tlat, tlon, tnorm, tyaw = getGTpose(C, qname)
    qzx = geom.lltom(olat,olon,tlat,tlon)

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

    # Solve for query pose using constrained homography or essential matrix
    runflag = param['runflag']
    if method == 'hom' and runflag != 6:
        runflag = 0
        parameters = ( Rq, Rd, cyaw, np.nan, runflag, param['inlier_error'], param['ransac_iter'] )
        matches, pose = solveGeom(matches,parameters)
    elif method == 'ess':
        runflag = 5
        parameters = ( Rq, Rd, cyaw, np.nan, runflag, param['inlier_error'], param['ransac_iter'] )
        matches, pose = solveGeom(matches,parameters)
    elif method == 'hom':
        pr = geom.YPRfromR(Rq)[1:]
        dyaws = range(-30,30,5)
        yaw_matches, yaw_poses, yerrs = [], [], []
        for dyaw in dyaws:
            yaw_Rq = geom.RfromYPR(cyaw+dyaw, pr[0], pr[1])
            parameters = ( yaw_Rq, Rd, cyaw+dyaw, np.nan, 2, param['inlier_error'], param['ransac_iter'] )
            yaw_matches.append(matches.copy())
            print 'Query yaw set to %d degrees' % (cyaw+dyaw)
            m, p = solveGeom(matches,parameters)
            yaw_poses.append(p)
            yerrs.append(m['ierr'])
        bestidx = np.argmin(yerrs)
        matches, pose = yaw_matches[bestidx], yaw_poses[bestidx]

    # extract pose parameters
    tray = pose[:3]
    qbear = pose[3]
    nbear = pose[4] if method=='hom' else np.nan
    

    # Get scaled translation for query location
    wRq_pr = geom.YPRfromR(Rq)[1:]
    comp_wRq = geom.RfromYPR(qbear, wRq_pr[0], wRq_pr[1])
    qloc = scalefrom3d(matches, tray, comp_wRq)[0]

    # temporarily get yaw error
    qyaw_error = int(abs(np.mod(tyaw-qbear,360)))
    qyaw_error = qyaw_error if qyaw_error<180 else abs(qyaw_error-360)

    # compute location errors wrt estimated query locations
    loc_err = ( (qloc[0]-qzx[1])**2 + (qloc[2]-qzx[0])**2 )**0.5
    gps_err = ( (gzx[1]-qzx[1])**2 + (gzx[0]-qzx[0])**2 )**0.5

    # compute the angle difference between T and ground truth translation
    tyaw_error = int(abs( 180/np.pi * np.arccos( np.abs(qloc[0]*qzx[1]+qloc[2]*qzx[0]) / (alg.norm([qloc[0],qloc[2]])*alg.norm(qzx)) ) ))

    # compute the plane normal angle error
    nyaw_error = np.nan if np.isnan(nbear) or np.isnan(tnorm) else np.mod(int(abs(nbear-tnorm)),180)
    nyaw_error = nyaw_error if nyaw_error<90 else abs(nyaw_error-180)

    # write pose estimation results to file
    yaw_err = np.nan
    with open(param['pose_file'],'a') as pose_file:
        print >>pose_file, '\t'.join([qname, str(loc_err), str(gps_err), \
            str(tyaw_error), str(qyaw_error), str(nyaw_error), str(matches['nvrf']), \
            str(matches['numq']), str(matches['nmat']), str(matches['ierr']), \
            str(qloc[0]), str(qloc[2]), str(norm_err), str(yaw_err), str(runflag)])

    # draw matches
    close = loc_err < 10
    imgpath = os.path.join( param['resultsdir'] , qname + ';gt' + str(gtStatus) + ';' + 'close' + str(close) + ';locerr=' + str(loc_err) + ';tAng=' + str(tyaw_error) + ';qAng=' + str(qyaw_error) + ';nAng=' + str(nyaw_error) + ';' + dbmatch + '.jpg')
    draw_matches(matches, qimg, dbimg, imgpath)

    print 'Computed yaw / ground truth yaw        : %d / %d' % (qbear,tyaw)
    if method == 'hom': print 'Computed normal bearing / ground truth : %d / %d' % (nbear,tnorm)
    print 'Computed query location relative to db     : %.1f, %.1f, %.1f' % tuple(qloc)
    print 'Ground truth query location relative to db : %.1f,  - , %.1f' % (qzx[1],qzx[0])

    input = (Rq, Rd)

    return qloc, loc_err, matches, input


def solveGeom(matches, parameters):
    Rq, Rd, qyaw, nyaw, runflag, inlerr, iter = parameters
    if runflag < 4: # homography
        matches, pose = solveHomography.constrainedHomography( matches , \
            Rd , Rq , qyaw , nyaw , runflag , maxerr=inlerr , maxiter=iter)
        if matches['nvrf']==0:
            matches, pose = solveHomography.constrainedHomography( matches , \
                Rd , Rq , qyaw , nyaw , runflag , maxerr=2*inlerr , maxiter=iter)
        if matches['nvrf']==0:
            matches, pose = solveHomography.constrainedHomography( matches , \
                Rd , Rq , qyaw , nyaw , runflag , maxerr=5*inlerr , maxiter=iter)
    else: # essential matrix
        matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
            Rd , Rq , qyaw , runflag , maxerr=inlerr , maxiter=iter)
        if matches['nvrf']==0:
            matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
                Rd , Rq , qyaw , runflag , maxerr=2*inlerr , maxiter=iter)
        if matches['nvrf']==0:
            matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
                Rd , Rq , qyaw , runflag , maxerr=5*inlerr , maxiter=iter)
    return matches, pose


def scalefrom3d(matches, tray, wRq):

    # extract inliers
    imask = matches['imask']
    qray = matches['qray'][imask,:]
    w3d = matches['w3d'][imask,:]

    # compute direction of 3d point from query image
    wray = tp(np.dot(wRq,tp(qray)))

    # set up Aq=k equation representing intersection of 2 lines
    range_numi = range(np.sum(imask))
    A = [ np.array([[ wray[i,2], -wray[i,0] ], [ tray[2], -tray[0] ]]) for i in range_numi ]
    k = [ np.array( [ wray[i,2]*w3d[i,0]-wray[i,0]*w3d[i,2] , 0 ] ) for i in range_numi ]

    # solve for intersection of 2 lines to get query location
    t_int = [ alg.solve(A[i],k[i]) for i in range_numi ]
    
    # compute the corresponding scale factors attached to y
    idx = int(tray[2]>tray[0])
    scales = [ t_int[i][idx]/tray[2*idx] for i in range_numi ]

    # find best collection of scale factors
    affprm = [ 1 , 0.1 ] # ransac loop chooses all scale factors within affprm[0]+s*affprm[1] meters of chose scale factor s
    bnum, bmask = 0, np.bool_(np.zeros(len(scales)))
    for i in range(len(scales)):
        s = scales[i]
        mask = np.abs(s-scales) < affprm[0] + affprm[1] * np.abs(s)
        if np.sum(mask) > bnum : bnum, bmask = np.sum(mask), mask
    tdist = np.mean(scales[bmask])
    t = tdist*tray

    return t, tdist