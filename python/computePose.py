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
#import vanPts

MAX_PIXEL_DEVIATION = 1
FALLBACK_PIXEL_DEVIATIONS = [2,1]
CONFIDENCE_LEVEL = .99999
BAD_HOMOGRAPHY_DET_THRESHOLD = .005
ROT_THRESHOLD_RADIANS = 0.2 # .1 ~ 5 deg


def highresSift(C, Q, dbsift):

    start = time.time()
    maxmatch, maxdist, maxangle = 5, 5000000, np.pi#/6
    idx = dbsift.rfind('/')
    hrRematchFile = dbsift[:idx] + '/highres/rematch/' + Q.name + ';' + \
                    dbsift[idx+1:-8] + ';maxmatch=' + str(maxmatch) + \
                    ';maxdist=' + str(maxdist/1000) + 'k;maxangle=' + \
                    str(int(round(180/np.pi*maxangle))) + '.npz'
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
        qsift = Q.siftpath
        idx = qsift.rfind('/')
        qsift = ''.join([qsift[:idx],'/hrSIFT',qsift[idx:]])
        q = rdr.load_file(qsift)
        idx = dbsift.rfind('/')
        dbsift = ''.join([dbsift[:idx],'/highres/hrSIFT',dbsift[idx:]])
        db = rdr.load_file(dbsift)
        flann = pyflann.FLANN()
        results, dists = flann.nn(db['vec'], q['vec'], maxmatch, algorithm='linear')
        if maxmatch == 1:
            results = [[r] for r in results]
            dists = [[d] for d in dists]
        matches = {'qidx': np.array([]), \
                   'q2d': np.zeros([0,2]), \
                   'qprm': np.zeros([0,2]), \
                   'd2d': np.zeros([0,2]), \
                   'dprm': np.zeros([0,2]), \
                   'nnd': np.array([]) }
        nmat, numq = 0, 0
        for i in xrange(len(results)):
            grads = np.mod(np.array([q[i]['geom'][3]-db[results[i][k]]['geom'][3] for k in range(maxmatch)]),2*np.pi)
            idx = np.nonzero( np.logical_and( np.array(dists[i])<maxdist , np.logical_or( grads<maxangle , grads-2*np.pi>-maxangle ) ) )[0]
            qadd = len(idx)
            if qadd == 0:
                continue
            q2d = np.tile(np.array([q[i]['geom'][1],q[i]['geom'][0]]),[qadd,1])
            qprm = np.tile(np.array(q[i]['geom'][2:].copy()),[qadd,1])
            d2d = np.array( [ [ db[results[i][k]]['geom'][1] , db[results[i][k]]['geom'][0] ] for k in idx ] )
            dprm = np.array( [ db[results[i][k]]['geom'][2:].copy() for k in idx ] )
            nnd = np.array( [ dists[i][k] for k in idx ] )
            matches['q2d'] = np.append(matches['q2d'],q2d,axis=0)
            matches['qprm'] = np.append(matches['qprm'],qprm,axis=0)
            matches['d2d'] = np.append(matches['d2d'],d2d,axis=0)
            matches['dprm'] = np.append(matches['dprm'],dprm,axis=0)
            matches['nnd'] = np.append(matches['nnd'],nnd,axis=0)
            matches['qidx'] = np.append(matches['qidx'],nmat)
            nmat += qadd
            numq += 1
        matches['qidx'] = np.append(matches['qidx'],nmat)
        matches['nmat'] = nmat
        matches['numq'] = numq
        np.savez( hrRematchFile, qidx=matches['qidx'], q2d=matches['q2d'], qprm=matches['qprm'],
                  d2d=matches['d2d'], dprm=matches['dprm'], nnd=matches['nnd'] )
        # END OF RUNNING NN SEARCH OR LOADING NN FILE
    # fill in other keys
    nmat = matches['nmat']
    numq = matches['numq']
    matches['qray'] = np.zeros([nmat,3])
    matches['dray'] = np.zeros([nmat,3])
    matches['w3d'] = np.zeros([nmat,3])
    matches['hmask'] = np.zeros(nmat)
    matches['herr'] = -1*np.ones(nmat)
    matches['hvrf'] = 0
    # Print rematch statistics
    print 'Number of query features matched: %d' % numq
    print 'Total number of feature matches: %d' % nmat
    print 'Average number of db matches considered = %.1f' % (float(nmat)/numq)
    print 'High res rematch took %.1f seconds.' % (time.time() - start)
    return matches


def lowresSift(C, Q, dbsift):
    
    start = time.time()
    rdr = reader.get_reader(C.params['descriptor'])
    qsift = Q.siftpath
    q = rdr.load_file(qsift)
    db = rdr.load_file(dbsift)
    flann = pyflann.FLANN()
    results, dists = flann.nn(db['vec'], q['vec'], 1, algorithm='linear')
    results, dists, q = np.array(results), np.array(dists), np.array(q)
    idx = np.argsort(dists)
    results = results[idx]
    dists = dists[idx]
    q = q[idx]
    count = 10000
    matches = []
    closed = set()
    for i in range(0,len(results)):
        if results[i] not in closed and dists[i] < 40000:
            closed.add(results[i])
            atom = {'db': db[results[i]]['geom'].copy(),
                      'query': q[i]['geom'].copy()}
            matches.append(atom)
            count -= 1
        if count == 0:
            break
    print 'Low res rematch took %.1f seconds.' % (time.time()-start)
    return matches


def draw_matches(C, Q, matches, db_img, out_img):

    print 'Drawing homography verified sift matches...'

    def scaledown(image, max_height):
        scale = 1.0
        hs = float(image.size[1]) / max_height
        if hs > 1:
            w,h = image.size[0]/hs, image.size[1]/hs
            scale /= hs
            image = image.resize((int(w), int(h)), Image.ANTIALIAS)
        return image, scale

    assert os.path.exists(db_img)
    a = Image.open(Q.jpgpath)
    b = Image.open(db_img)
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

    match_idx = np.nonzero( np.ones(matches['nmat'] ) )[0] # np.nonzero(matches['hmask'])[0]
    for idx in match_idx:
        start = scale*matches['q2d'][idx,:]
        stop = matches['d2d'][idx,:]
        stop[0] += off
        try:
            xdrawcircle(start,'red')
        except ValueError:
            print start
#        xdrawline((start,stop),'orange')
        xdrawcircle(stop,'red')

    target.save(out_img, 'jpeg', quality=90)


def createHR3dMap(hr3dmap, map3d, dbs, size):
    map3dhr = np.zeros([size[0],size[1],3])
    mapkeys = map3d.keys()
    infarr3 = arr([np.inf,np.inf,np.inf])
    for i in xrange(size[0]):
        print 'Processing row %d of 3d map...' % i
        for j in xrange(size[1]):
            d3d = map3d[tuple(mapkeys[np.argmin([alg.norm(dbs*arr([i,j])-k) for k in mapkeys])])]
            if d3d is None:
                map3dhr[i,j,:] = infarr3
            else:
                map3dhr[i,j,:] = [d3d['lat'],d3d['lon'],d3d['alt']]
    np.save(hr3dmap,map3dhr)
    return map3dhr


def pose_triangle(C, Q, dbimg, dbsift, udir, vpN=False):

    # get high res sift and img paths
    idx = dbimg.rfind('/')
    hrdbimg = dbimg[:idx] + '/highres' + dbimg[idx:]
    path3dmap = dbimg[:idx] + '/highres/hr3dmap' + dbimg[idx:-4] + '.map'

    # get high res sift rematch
    matches = highresSift(C, Q, dbsift)

#    # Get image information
#    qsource = render_tags.QueryImageInfo(Q.datasource)
#    info = os.path.join(C.infodir, os.path.basename(hrdbimg)[:-4] + '.info')
#    dbsource = render_tags.EarthmineImageInfo(hrdbimg, info)
#    lrinfo = os.path.join(C.infodir, os.path.basename(dbimg)[:-4] + '.info')
#    lrdb = render_tags.EarthmineImageInfo(dbimg, lrinfo)
#
#    # Set Kq, Rq
#    wx,wy = qsource.image.size
#    fov = qsource.view_angle[0]
#    Kq = geom.cameramat(wx, wy, fov)
#    Kqinv = alg.inv(Kq)
#    y,p,r = qsource.yaw, qsource.pitch, qsource.roll
#    Rq = geom.RfromYPR(y,p,r) # camera orientation (camera to world)
#
#    # Set Kd, Rd
#    wx,wy = dbsource.image.size
#    fov = dbsource.fov
#    Kd = geom.cameramat(wx, wy, fov)
#    Kdinv = alg.inv(Kd)
#    lr_matches = lowresSift(C, Q, dbsift)
#    Rd,foo = pnp.dbsolve(C, Rq, lr_matches, dbsift, dbimg)
#    db_err = alg.norm(np.array(foo))
#
#    # Fill out match information
#    print 'Filling out query specific information on sift matches...'
#    start = time.time()
#    nmat = matches['nmat']
#    matches['qray'] = tp(dot(Kqinv,np.append(tp(matches['q2d']),[np.ones(nmat)],0)))
#    matches['dray'] = tp(dot(Kdinv,np.append(tp(matches['d2d']),[np.ones(nmat)],0)))
#    lr3dmap = C.pixelmap.open(dbsift)
#    lrKeys = lr3dmap.keys()
#    if os.path.isfile(path3dmap):
#        saveMap = False
#        file3dmap = open(path3dmap,'r')
#        try:
#            hr3dmap = pickle.load(file3dmap)
#        except EOFError:
#            saveMap = True
#            hr3dmap = lr3dmap
#        file3dmap.close()
#    else:
#        saveMap = True
#        hr3dmap = lr3dmap
#    dbs = float(lrdb.image.size[0]) / dbsource.image.size[0]
#    olat,olon,oalt = dbsource.lat,dbsource.lon,dbsource.alt # database location
#    iremove = np.array([])
#    for i in xrange(nmat):
#        key = (int(round(dbs*matches['d2d'][i,0])),int(round(dbs*matches['d2d'][i,1])))
#        if key in hr3dmap:
#            data3d = hr3dmap[key]
#        else:
#            saveMap = True
#            data3d = lr3dmap[tuple(lrKeys[np.argmin([alg.norm(arr(key)-k) for k in lrKeys])])]
#            hr3dmap[key] = data3d
#        if data3d is None:
#            iremove = np.append(iremove,i)
#        else:
#            zx = geom.lltom(olat,olon,data3d['lat'],data3d['lon'])
#            matches['w3d'][i,:] = np.array([zx[1],oalt-data3d['alt'],zx[0]])
#    # Remove indices without 3d information
#    matches['q2d'] = np.delete(matches['q2d'],iremove,0)
#    matches['qprm'] = np.delete(matches['qprm'],iremove,0)
#    matches['qray'] = np.delete(matches['qray'],iremove,0)
#    matches['d2d'] = np.delete(matches['d2d'],iremove,0)
#    matches['dprm'] = np.delete(matches['dprm'],iremove,0)
#    matches['dray'] = np.delete(matches['dray'],iremove,0)
#    matches['w3d'] = np.delete(matches['w3d'],iremove,0)
#    matches['nnd'] = np.delete(matches['nnd'],iremove,0)
#    matches['nmat'] -= len(iremove)
#    matches['qidx'] = np.array(list(set( [ qi-np.sum(iremove<qi) for qi in matches['qidx']] )))
#    matches['numq'] = len(matches['qidx'])-1
#    matches['hmask'] = np.zeros(matches['nmat'])
#    matches['herr'] = -1 * np.ones(matches['nmat'])
#    matches['hvrf'] = 0
#    # save hr3dmap if changed
#    if saveMap:
#        file3dmap = open(path3dmap,'w')
#        pickle.dump(hr3dmap,file3dmap)
#        file3dmap.close()
#    print 'Retrieving 3d points took %.1f seconds.' % (time.time()-start)
#
#    # Get estimated groundd truth query location and normal direction
#    qlocs = open('/media/DATAPART2/ah/query5horizontal-locs.txt','r')
#    qstr = qlocs.read()
#    qlocs.close()
#    qstr = qstr[qstr.find(qsource.name):]
#    qstr = qstr[qstr.find('\t')+1:]
#    qidx = qstr.find('\t')
#    qlat = float(qstr[:qidx])
#    qstr = qstr[qidx+1:]
#    qidx = qstr.find('\t')
#    qlon = float(qstr[:qidx])
#    qstr = qstr[qidx+1:]
#    qidx = qstr.find('\n')
#    qbear = int(qstr[:qidx])
#    qzx = geom.lltom(olat,olon,qlat,qlon)
#
#    # Retrieve GPS query location
#    glat,glon = qsource.lat,qsource.lon
#    gzx = geom.lltom(olat,olon,glat,glon)

#    # Solve for normal vector using vanishing points
#    vpN = True
#    if vpN:
#        vbear, vbErr = computeNorm(Q.jpgpath,hrdbimg,Rq,Rd,Kq,Kd)
#    else:
#        vbear = np.nan
#    vbeardeg = 180/np.pi*vbear - 180*round((vbear*180/np.pi-qbear)/180.)
#    vAng_err = np.nan if np.isnan(vbear) else int(abs(vbeardeg-qbear))
#    print vAng_err

#    # Solve for query pose using constrained homography
#    mxit = 10000
#    matches, ct, nbear = constrainedHomography(matches,Rd,Rq,nbear=vbear,maxerr=.03,maxiter=mxit)
#    if matches['hvrf']==0 or alg.norm(ct)>50:
#        print 'Attempting to solve constrained homography with relaxed parameters.'
#        matches, ct, nbear = constrainedHomography(matches,Rd,Rq,nbear=vbear,maxerr=.06,maxiter=mxit)
#    if matches['hvrf']==0 or alg.norm(ct)>50:
#        print 'Attempting to solve constrained homography with even more relaxed parameters.'
#        matches, ct, nbear = constrainedHomography(matches,Rd,Rq,nbear=vbear,maxerr=.12,maxiter=mxit)
#    nbear = np.mod( 180 + 180/np.pi * nbear , 360 )
#
#    # compute location errors wrt estimated query locations
#    ch_err = ( (ct[0]-qzx[1])**2 + (ct[2]-qzx[0])**2 )**0.5
#    gps_err = ( (gzx[1]-qzx[1])**2 + (gzx[0]-qzx[0])**2 )**0.5
#
#    # compute the angle difference between T and ground truth translation
#    tAng_err = int(abs( 180/np.pi * np.arccos( (ct[0]*qzx[1]+ct[2]*qzx[0]) / (alg.norm([ct[0],ct[2]])*alg.norm(qzx)) ) ))
#
#    # compute the plane normal angle error
#    nAng_err = np.nan if np.isnan(nbear) else int(abs(nbear-qbear))
#    vAng_err = np.nan if np.isnan(vbear) else int(abs(vbear-qbear))
#
#    # write pose estimation results to file
#    with open(C.reproj_file,'a') as rp_file:
#        print >>rp_file, '\t'.join([os.path.split(Q.jpgpath)[1][:-4], str(ch_err), str(gps_err),
#              str(tAng_err), str(nAng_err), str(matches['hvrf']), str(matches['nmat']), str(ct[0]), str(ct[2]), str(db_err)])
#
    # draw matches
#    close = ch_err < 10
#    imgpath = udir + '/' + Q.name + ';gtTrue;' + 'close' + str(close) + ';err=' + str(ch_err) + ';tAng=' + str(tAng_err) + ';nAng=' + str(nAng_err) + ';feature_pairs.jpg'
    imgpath = udir + '/' + Q.name + '.jpg'
    draw_matches(C, Q, matches, hrdbimg, imgpath)

#    # write pose estimation results to file
#    with open(C.reproj_file,'a') as rp_file:
#        print >>rp_file, '\t'.join([os.path.split(Q.jpgpath)[1][:-4], str(vAng_err)])


    return ct, ch_err, matches


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
    knownN = np.isnan(nbear)

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
            verr = np.reshape(homerrf(mp,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2])>stoperr).any():
                mp[:3] = dot(wRd,-k*r)
                verr = np.reshape(homerrf(mp,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2])>stoperr).any():
                continue
            mp = opt.leastsq(homerrf,mp,args=([mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),warning=False)[0]
            verr = np.reshape(homerrf(mp,[mq1,mq2],[md1,md2],[mw1,mw2],dRq,wRd,spd),[-1,3])
            if (abs(verr[:,:2]>stoperr)).any() or (abs(verr[:,2])>spd).any(): # verify homography solution and 3d world points
                continue
            errs = np.sum(np.reshape(hom_errf(mp,matches['qray'],matches['dray'],matches['w3d'],dRq,wRd,spd),[-1,3])**2,1)**0.5
        emask = errs < maxerr
        imask = np.bool_(np.zeros(matches['nmat']))
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
            ip = opt.leastsq(homerrf,mp,args=(iq,id,iw,dRq,wRd,spd),warning=False)[0]
            ierr = alg.norm(homerrf(ip,iq,id,iw,dRq,wRd,spd)) / numi
            if ierr < berr:
                berr = ierr
                T = ip[4]*ip[:3]
                nbear = ip[3]
                matches['Hprm'] = ip
                matches['Hmat'] = dot(tp(dRq),np.eye(3,3)-dot(r,tp(n)))
                matches['Tvec'] = T
                matches['Pnorm'] = arr([np.sin(ip[3]),0,np.cos(ip[3])])
                matches['hmask'] = imask
                matches['herr'][imask,:] = np.sum(np.reshape(homerrf(ip,iq,id,iw,dRq,wRd,spd),[-1,3])**2,1)**0.5
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
    err = arr( [ [ q[i,0]-Hd[i,0]/Hd[i,2] , q[i,1]-Hd[i,1]/Hd[i,2] , spd*(prm[3]-dot(wn,w[0,:])) ] for i in range(len(q)) ] )
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
    k = m[0] / np.cos(e1-f[0])
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

#def essmat_fromR(yd, yq, dR):
#
#    ncand = len(yd)
#    if ncand < 3:
#        raise Exception('Too few correspondences to compute essential matrix.')
#    x, y = arr(yd), arr(yq)
#    z = arr([dot(tp(dR),q) for q in y])
#
#    # RANSAC outer loop
#    infarr = arr([np.inf,np.inf,np.inf])
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
#            M[i,:] = arr( [ mx[i,1]*mz[i,2]-mx[i,2]*mz[i,1] ,
#                            mx[i,2]*mz[i,0]-mx[i,0]*mz[i,2] ,
#                            mx[i,0]*mz[i,1]-mx[i,1]*mz[i,0] ] )
#        u,s,v = alg.svd(M)
#        mt = v[2,:]
#        mT = arr([[0,-mt[2],mt[1]],
#                  [mt[2],0,-mt[0]],
#                  [-mt[1],mt[0],0]])
#
#        err = arr([dot(z[i,:],dot(mT,x[i,:])) / (alg.norm(z[i,:])*alg.norm(x[i,:])) for i in range(ncand)])
#        inmask = abs(err) < maxerr
#        ix, iz, im = np.compress(inmask,x,0), np.compress(inmask,z,0), sum(inmask)
#
#        if sum(inmask) >= minfit:
#            M = np.zeros([im,3])
#            for i in range(im):
#                M[i,:] = arr( [ ix[i,1]*iz[i,2]-ix[i,2]*iz[i,1] ,
#                                ix[i,2]*iz[i,0]-ix[i,0]*iz[i,2] ,
#                                ix[i,0]*iz[i,1]-ix[i,1]*iz[i,0] ] )
#            u,s,v = alg.svd(M)
#            ierr = (s[2]/s[0])**0.5 / im
#            if ierr < berr:
#                bmask, berr = inmask, ierr
#                bt = v[2,:]
#                bT = arr([[0,-bt[2],bt[1]],
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
#    pnorm = arr([p/alg.norm(p) for p in pnorm]) #normalize
##    print 'Plane normals'
##    print pnorm
#
#    # RANSAC outer loop
#    iter, maxiter = 0, 1000
#    bt, bmask, berr = arr([np.inf,np.inf,np.inf]), np.ones(npts), np.inf
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
#    A = arr([[y[1],-y[0]],[z[1],-z[0]]])
#    k = y[1]*x[0]-y[0]*x[1]
#    k = arr([k,0])
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
#    q2d = [arr([m[0],m[1]]) for m in matches]
#    db2d = [arr([m[2],m[3]]) for m in matches]
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
#    #print arr([d['alt'] for d in db3d])
#    #print dbsource.alt
#    y = [dbsource.alt-d['alt'] for d in db3d]
#    xyz = [[zx[i][1],y[i],zx[i][0]] for i in range(len(y))]
#    #print -tp(arr(y))
#
#    # Set K, e3, Rhat and get list of ci, xi
#    wx,wy = qsource.pgmsize[0], qsource.pgmsize[1]
#    tx,ty = qsource.view_angle[0], qsource.view_angle[1]
#    imsize = arr([wx,wy])
#    f1, f2 = (wx-1)/(2*np.tan(tx/2)), (wy-1)/(2*np.tan(ty/2))
#    f = (f1+f2) / 2
#    print (wx,wy,tx,f)
#    lsprm = np.concatenate((arr([f]),imsize),1)
#    K = arr([[f,0,(wx-1)/2.0],
#             [0,f,(wy-1)/2.0],
#             [0,0,1]])
#    e3 = tp(arr([[0,0,1]]))
#    y,p,r = qsource.yaw, qsource.pitch, qsource.roll
#    #r = -np.pi/2
#    print [180*y/np.pi,180*p/np.pi,180*r/np.pi]
#    Ry = arr([[np.cos(y),0,-np.sin(y)],
#              [0,1,0],
#              [np.sin(y),0,np.cos(y)]])
#    Rx = arr([[1,0,0],
#              [0,np.cos(p),np.sin(p)],
#              [0,-np.sin(p),np.cos(p)]])
#    Rz = arr([[np.cos(r),np.sin(r),0],
#              [-np.sin(r),np.cos(r),0],
#              [0,0,1]])
#    Rhat = dot(Ry,dot(Rx,Rz))
#    ci = [tp(arr([[q[0],q[1],1]])) for q in q2d]
#    xi = [tp(arr([[x[0],x[1],x[2]]])) for x in xyz]
#
#    # Minimize linearized error for T based on orientation R from sensors
#    # Note that T represents translation from camera to global
#    # Global to camera is t = -inv(R)T or T = -Rt
#    # Solution to minimization of linear error is T = -inv(A)*x
#    N = len(xi)
#    A = arr([[0,0,0],
#             [0,0,0],
#             [0,0,0]])
#    x = tp(arr([[0,0,0]]))
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
#    ci = tp(np.reshape(arr(ci)[:,:2,0],(numpts,-1)))
#    xi = tp(np.reshape(arr(xi),(numpts,-1)))
#
#    # init orientation
#    rot = arr([y,p,r])
#
#    # full optimization with least squares
#    errfunc = reproj_errfunc1
#    p0 = arr([y,p,r,T[0],T[1],T[2]])
#    p0 = arr([0,0,0])
#    p2, success = opt.leastsq(errfunc,p0,args=(rot,xi,ci,lsprm))
#
#    # optimization of T only with ransac least squares
#    errfunc = reproj_errfunc1
#    p0 = arr([0,0,0])
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
#    #T = tp(arr(p[0:3]))
#    ry,rx,rz = p[0:3]
#    T = tp(arr(p[3:6]))
#    Ry = arr([[np.cos(ry),0,-np.sin(ry)],
#              [0,1,0],
#              [np.sin(ry),0,np.cos(ry)]])
#    Rx = arr([[1,0,0],
#              [0,np.cos(rx),np.sin(rx)],
#              [0,-np.sin(rx),np.cos(rx)]])
#    Rz = arr([[np.cos(rz),np.sin(rz),0],
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
#    K = arr([[f,0,(imsize[0]-1)/2.0],
#             [0,f,(imsize[1]-1)/2.0],
#             [0,0,1]])
#    Ry = arr([[np.cos(ry),0,-np.sin(ry)],
#              [0,1,0],
#              [np.sin(ry),0,np.cos(ry)]])
#    Rx = arr([[1,0,0],
#              [0,np.cos(rx),np.sin(rx)],
#              [0,-np.sin(rx),np.cos(rx)]])
#    Rz = arr([[np.cos(rz),np.sin(rz),0],
#              [-np.sin(rz),np.cos(rz),0],
#              [0,0,1]])
#    R = dot(Ry,dot(Rx,Rz))
#    tmp = np.concatenate((R,tp(arr([T]))),1)
#    MAT = dot(K,np.concatenate((R,tp(arr([T]))),1))
#    U,S,V = alg.svd(MAT)
#    pt = V[3][:]
#    print [pt[0]/pt[3], p[1]/pt[3], pt[2]/pt[3]]
#
#    # test some 3d to 2d projections
#    param = np.concatenate((p0,arr([f]),imsize),1)
#    pi = proj3dto2d(param,xi)
#    map = arr([[ci[0,i],ci[1,i],pi[0,i],pi[1,i]] for i in range(numpts)])
#    err = np.sqrt( (map[:,0]-map[:,2])**2 + (map[:,1]-map[:,3])**2 )
#    avg = np.average(err,0)
#    print 'Average pixel reprojection error for p0: {0}'.format(avg)
#
#    # test some 3d to 2d projections
#    param = np.concatenate((p,arr([f]),imsize),1)
#    pi = proj3dto2d(param,rscX)
#    map = arr([[rscC[0,i],rscC[1,i],pi[0,i],pi[1,i]] for i in range(np.shape(rscC)[1])])
#    err = np.sqrt( (map[:,0]-map[:,2])**2 + (map[:,1]-map[:,3])**2 )
#    avg = np.average(err,0)
#    print 'Average pixel reprojection error for inliers: {0}'.format(avg)
#
#    # test some 3d to 2d projections
#    param = np.concatenate((p,arr([f]),imsize),1)
#    pAll = proj3dto2d(param,xi)
#    mapAll = arr([[ci[0,i],ci[1,i],pAll[0,i],pAll[1,i]] for i in range(np.shape(ci)[1])])
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
#    T = tp(arr(param[3:6]))
#    f = param[6]
#    imsize = param[7:]
#    K = arr([[f,0,(imsize[0]-1)/2.0],
#             [0,f,(imsize[1]-1)/2.0],
#             [0,0,1]])
#    Ry = arr([[np.cos(ry),0,-np.sin(ry)],
#              [0,1,0],
#              [np.sin(ry),0,np.cos(ry)]])
#    Rx = arr([[1,0,0],
#              [0,np.cos(rx),np.sin(rx)],
#              [0,-np.sin(rx),np.cos(rx)]])
#    Rz = arr([[np.cos(rz),np.sin(rz),0],
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
