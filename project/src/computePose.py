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
import numpy.random as rnd
from numpy import transpose as tp
import geom
import os
import solveHomography
import solveEssMatrix
import vp_analysis
import pickle

MAX_PIXEL_DEVIATION = 1
FALLBACK_PIXEL_DEVIATIONS = [2,1]
CONFIDENCE_LEVEL = .99999
BAD_HOMOGRAPHY_DET_THRESHOLD = .005
ROT_THRESHOLD_RADIANS = 0.2 # .1 ~ 5 deg


def highresSift(C, Q, dbmatch):

    # timing info
    start = time.time()

    # set sift paths
    qname = Q.name
    dbsift = os.path.join(C.hiresdir,dbmatch+'sift.txt')
    qsift = os.path.join(C.querydir,'hires',qname+'sift.txt')

    # high res rematch parameters
    maxmatch, maxangle, maxdist = 1, np.pi/3, 10**7
    maxratio = C.pose_param['maxratio']

    # rematch file
    filename = qname + ';' + dbmatch + ';maxratio=' + str(maxratio) + \
                    ';maxmatch=' + str(maxmatch) + ';maxdist=' + \
                    str(maxdist/1000) + 'k;maxangle=' + \
                    str(int(round(180/np.pi*maxangle))) + '.npz'
    hrRematchFile = os.path.join( C.querydir, 'hires', 'siftmatch', filename )

    ### HIGH RES REMATCH ###
    matches = {}
    if not os.path.isdir(os.path.dirname(hrRematchFile)): os.path.mkdir(os.path.dirname(hrRematchFile))
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
        results, dists = flann.nn(db['vec'], q['vec'], 1+maxmatch, algorithm='linear')
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
        print 'High res rematch took %.1f seconds.' % (time.time() - start)
        np.savez( hrRematchFile, qidx=matches['qidx'], d2d=matches['d2d'], dprm=matches['dprm'],
                  q2d=matches['q2d'], qprm=matches['qprm'], nnd=matches['nnd'] )
        # END OF RUNNING NN SEARCH OR LOADING NN FILE
    # fill in other keys
    nmat = matches['nmat']
    numq = matches['numq']
    matches['dray'] = np.zeros([nmat,3])
    matches['qray'] = np.zeros([nmat,3])
    matches['w3d'] = np.zeros([nmat,3])
    matches['ddep'] = np.zeros(nmat)
    matches['plane'] = np.int_( -1 * np.ones(nmat) )
    matches['weight'] = np.zeros(nmat)
    matches['imask'] = np.bool_(np.zeros(nmat))
    matches['hvrf'] = 0
    # Print rematch statistics
    print 'Number of query features matched: %.0f' % numq
    print 'Total number of feature matches: %.0f' % nmat
    print 'Average number of database matches considered = %.1f' % (float(nmat)/numq)
    return matches


def draw_dbimage(C, Q, matchedimg, match):

    print 'Drawing query and database images side by side...'

    qname = Q.name
    dbimg = os.path.join(C.hiresdir,matchedimg+'.jpg')
    qimg = os.path.join(C.querydir,'hires',qname+'.jpg')
    outimg = os.path.join(C.pose_param['resultsdir'],qname+';'+str(match)+';'+matchedimg+'.jpg')

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

    target.save(outimg, 'jpeg', quality=90)


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

    match_idx = np.nonzero(matches['imask'])[0]
    #match_idx = np.nonzero( np.zeros(matches['nmat'] ) )[0]
    print 'Number of inliers / total correspondences : %d / %d' % (matches['numi'],matches['nmat'])
    for idx in match_idx:
        start = scale*matches['q2d'][idx,:]
        stop = matches['d2d'][idx,:]
        stop[0] += off
        xdrawcircle(start,'red')
        xdrawline((start,stop),'orange')
        xdrawcircle(stop,'red')

    target.save(outimg, 'jpeg', quality=90)


def match_info(C, Q, matches, dbmatch, wRd, domplane=-1):

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
#        depthfail = depth == 0 or depth == 65535
        depthfail = depth == 0 or depth > 10000 # 100+ meters away eliminated
        planefail = domplane != -1 and plane != domplane
        groundfail = C.pose_param['remove_ground'] and plane == 0
        if depthfail or planefail or groundfail:
            iremove = np.append(iremove,i)
            continue
        matches['plane'][i] = plane if plane != 255 else -1
        depth = depth / 100.
        matches['ddep'][i] = depth
        matches['weight'][i] = 20/depth if ( depth>20 and C.pose_param['use_weight'] ) else 1
        wray = np.dot( wRd , matches['dray'][i,:] )
        matches['w3d'][i,:] = depth * wray
    # Remove indices without 3d information
    matches['domplane'] = domplane
    matches['q2d'] = np.delete(matches['q2d'],iremove,0)
    matches['qprm'] = np.delete(matches['qprm'],iremove,0)
    matches['qray'] = np.delete(matches['qray'],iremove,0)
    matches['d2d'] = np.delete(matches['d2d'],iremove,0)
    matches['dprm'] = np.delete(matches['dprm'],iremove,0)
    matches['dray'] = np.delete(matches['dray'],iremove,0)
    matches['w3d'] = np.delete(matches['w3d'],iremove,0)
    matches['ddep'] = np.delete(matches['ddep'],iremove,0)
    matches['nnd'] = np.delete(matches['nnd'],iremove,0)
    matches['nmat'] -= len(iremove)
    matches['qidx'] = np.sort(np.array(list(set( [ qi-np.sum(iremove<qi) for qi in matches['qidx']] ))))
    matches['numq'] = len(matches['qidx'])-1
    matches['imask'] = np.bool_(np.zeros(matches['nmat']))
    print 'Removed %.0f points, leaving %.0f feature matches' % (len(iremove), matches['nmat'])
    print 'Retrieving 3d points and planes took %.1f seconds.' % (time.time()-start)
    return matches


def getGTpose(C, Q):
    qname = Q.name
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


def find_dbplanes(C, Q, dbmatch, Kdinv, wRd):
    print 'Finding database planes...'
    if C.QUERY != 'oakland1': return ( np.zeros(0), np.zeros(0), np.zeros((0,5)) )
    planefile = os.path.join(C.hiresdir,'planes',dbmatch+'-planes.npz')
    if not os.path.isdir(os.path.dirname(planefile)): os.path.mkdir(os.path.dirname(planefile))
    if os.path.isfile(planefile):
        pf_info = np.load(planefile);
        planes, sizes, params = pf_info['planes'], pf_info['sizes'], pf_info['params']
        print 'Database plane information loaded from cache.'
    else:
        planedata = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-planes.png') ) )
        h,w = planedata.shape
        maxplane = np.max(planedata[planedata!=255])
        plane_count = np.array( [ float(np.sum(np.reshape(planedata==k+1,-1)))/(h*w) for k in range(maxplane) ] )
        planes = 1 + np.nonzero(plane_count>0.1)[0]
        sizes = plane_count[planes-1]
        if len(planes) == 0: return planes, sizes, np.zeros((0,5),np.float)
        print 'Computing database plane normals...'
        start = time.time()
        params = np.array( [ planefrom3d(C, Q, dbmatch, planes[i], Kdinv, wRd) for i in range(len(planes)) ] )
        mask = np.logical_not(np.isnan(params[:,-1]))
        planes, sizes, params = planes[mask], sizes[mask], params[mask,:]
        np.savez(planefile,planes=planes,sizes=sizes,params=params)
        print 'Computing plane normals took %.1f seconds.' % (time.time()-start)
        np.savez(planefile,planes=planes,sizes=sizes,params=params)
    return planes, sizes, params


def planefrom3d(C, Q, dbmatch, domplane, Kdinv, wRd):

    if domplane == -1: return np.nan * np.zeros(5)

    # get 3d points on plane
    planes = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-planes.png') ) )
    depths = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-depth.png') ) )
    y, x = np.nonzero(planes==domplane)
    npts = len(x)
    pray = geom.normalrows( tp( np.dot( wRd , np.dot( Kdinv , np.concatenate( ([x],[y],[np.ones(npts)]) , 0 ) ) ) ) )
    pdep = depths[y,x]/100.0
    p3d = np.append( geom.vecmul(pray,pdep) , tp(np.array([np.ones(len(pray))])) , 1 )
    xz_pts = p3d[:,[0,2,3]]

    # RANSAC solve
    threshold, g = 2, np.array([0,1,0]) # meters
    bprm, bnumi, bmask = np.zeros(3), 0, np.bool_(np.zeros(npts))
    for i in range(1000):
        i1 = rnd.randint(0,npts)
        i2 = rnd.randint(0,npts-1)
        i2 = i2 if i2<i1 else i2+1
        i3 = rnd.randint(0,npts-2)
        i3 = i3 if i3<min(i1,i2) else ( i3+1 if i3+1<max(i1,i2) else i3+2 )
        inlpts = xz_pts[[i1,i2,i3],:]
        prm = geom.smallestSingVector(inlpts)
        prm = prm / geom.vecnorm(prm[:2])
        prm = -prm if prm[2]<0 else prm
        errs = np.abs(np.inner(xz_pts,prm))
        inlmask = errs < threshold
        numi = np.sum(inlmask)
        if numi > bnumi and float(numi)/npts > 0.5: bprm, bmask, bnumi = prm, inlmask, numi
    prm, numi, mask = bprm, bnumi, bmask

    # guided matching
    for i in range(10):
        if numi == 0: break
        prm = geom.smallestSingVector(xz_pts[mask,:])
        prm = prm / geom.vecnorm(prm[:2])
        prm = -prm if prm[2]<0 else prm
        errs = np.abs(np.inner(xz_pts,prm))
        mask = errs < threshold
        numi = np.sum(mask)

    # get error
    err = np.mean(np.abs(np.inner(xz_pts[mask,:],prm)))

    return np.array([prm[0],0,prm[1],prm[2],err])


def combine_planes(runflag,vnorms,dplanes,dsizes,dprms):
    # make a list of vp planes and db planes
    pyaws, planes, pconfs = np.zeros(0), np.zeros(0,np.int), np.zeros(0)
    while len(vnorms) > 0 and runflag > 6:
        vyaw = vnorms[0]
        pyaws, pconfs = np.append(pyaws,vyaw), np.append(pconfs,0)
        planes = np.append(planes,-1)
        vnorms  = np.delete(vnorms,0)
    while len(dplanes) > 0:
        dyaw = np.mod(180/np.pi*np.arctan2(dprms[0,0],dprms[0,2]),360)
        if dprms[0,-1] < 3:
            pyaws, planes, pconfs = np.append(pyaws,dyaw), np.append(planes,dplanes[0]), np.append(pconfs,0)
        dplanes, dsizes, dprms = np.delete(dplanes,0), np.delete(dsizes,0), np.delete(dprms,0,0)
    return pyaws, planes, pconfs


def estimate_pose(C, Q, dbmatch, gtStatus=None):

    # settings
    param = C.pose_param
    runflag = param['runflag']
    np.seterr(all='ignore')
    Q.datafile = os.path.join(C.pose_param['resultsdir'],'data_'+Q.name+'.txt')
    open(Q.datafile,'w').close()

    #####-----    PRINT RUN DETAILS    -----#####
    run_info = os.path.join(param['resultsdir'],param['run_info'])
    open(run_info,'w').close()
    with open(run_info,'a') as ri:
        if runflag == 11:   print >>ri, 'Method: Yaw, planes from VPs. Scale computed with homography.'
        elif runflag == 10: print >>ri, 'Method: Yaw, planes from VPs. Scale computed after homography.'
        if param['cheat']:  print >>ri, 'Ground truth yaw and plane used (cheating).'
        print >>ri, 'Inlier base error threshold: %.3f' % param['inlier_error']
        print >>ri, 'Base iteration scale: %d' % param['ransac_iter']
    #####-----    PRINT RUN DETAILS    -----#####

    # get high res db image and sift paths
    dbinfo = os.path.join(C.hiresdir, dbmatch + '.info')
    dbimg = os.path.join(C.hiresdir,dbmatch+'.jpg')
    dbsift = os.path.join(C.hiresdir,dbmatch+'sift.txt')
    dbsource = render_tags.EarthmineImageInfo(dbimg, dbinfo)

    # Set Kd, wRd, and db position
    wx,wy = dbsource.image.size
    fov = dbsource.fov
    Kd = geom.cameramat(wx, wy, fov)
    Kdinv = alg.inv(Kd)
    y,p,r = dbsource.yaw, dbsource.pitch, dbsource.roll
    wRd = geom.RfromYPR(y,p,r) # db camera orientation (camera to world)
    olat,olon,oalt = dbsource.lat,dbsource.lon,dbsource.alt # database location

    # get high res query information
    qname = Q.name
    qimg = os.path.join(C.querydir,'hires',qname+'.jpg')
    qsift = os.path.join(C.querydir,'hires',qname+'sift.txt')
    qsource = render_tags.QueryImageInfo(Q.datasource)
    glat,glon = qsource.lat,qsource.lon
    gzx = geom.lltom(olat,olon,glat,glon)
    timename = qname[-12:-10]+qname[-9:-7]+qname[-6:-4]#+qname[-3:]

    # Set Kq
    wx,wy = qsource.image.size
    fov = qsource.fov
    Kq = geom.cameramat(wx, wy, fov)
    Kqinv = alg.inv(Kq)
    cyaw,p,r = qsource.yaw, qsource.pitch, qsource.roll # cyaw - cell phone yaw

    # get high res sift rematch
    matches = highresSift(C, Q, dbmatch)
    with open(Q.datafile,'a') as df:
        print >>df, 'Number of matches | number of queries | ratio: %.0f | %.0f | %.2f' % (matches['nmat'], matches['numq'], float(matches['nmat'])/matches['numq'])
        print >>df, ''
        
    # Get estimated ground truth query location and normal direction
    tlat, tlon, tnorm, tyaw = getGTpose(C, Q)
    qzx = geom.lltom(olat,olon,tlat,tlon)

    # get query yaw and plane yaw from vanishing point anaylsis
    yawforvp = tyaw if param['cheat'] else np.nan
    vyaw, vnorms = vp_analysis.getQNyaws(C, Q, qimg, dbimg, qsource, yawforvp) 

    # get dominant planes
    dplanes, psizes, planeprms = find_dbplanes(C, Q, dbmatch, Kdinv, wRd)
    
    # match vanishing point planes to database planes
    pyaws, planes, pconfs = combine_planes(runflag,vnorms,dplanes,psizes,planeprms)

    print 'VP and DB Planes: ' + str(np.int_(pyaws)) + ', ' + str(planes)

    with open(Q.datafile,'a') as df:
#        print >>df, 'Planes detected with vanishing points:'
        for i in range(len(pconfs)):
            perr = np.round(np.mod(pyaws[i]-tnorm,360))
            perr = perr if perr<180 else 360-perr
            print >>df, 'Plane Yaw | DB plane | Confidence | Error : %3.0f | %d | %.2f | %.0f' % (pyaws[i],0 if planes[i]<0 else planes[i],pconfs[i],perr)
        yerr = np.round(np.mod(vyaw-tyaw,360))
        yerr = yerr if yerr<180 else 360-yerr
        print >>df, 'VP Yaw | Confidence | Error : %3.0f | %.2f | %.0f' % (vyaw,np.nan,yerr)
        print >>df, 'Cell yaw | True yaw | Plane : %3.0f | %3.0f  | %3.0f' % (cyaw,tyaw,tnorm)
        print >>df, ''

    # Set yaw value to be used
    if runflag >= 10: # vanishing point methods
        if np.isnan(vyaw): yaw, yawerr = cyaw, 0
        else: yaw, yawerr = vyaw, 0
    else: yaw, yawerr = cyaw, 0 # set cell phone yaw to use, plane normal
    wRq = geom.RfromYPR(yaw,p,r) # camera orientation (camera to world)

    ### --- THIS IS FOR CHEATING --- ###
    if param['cheat']:
        if not np.isnan(tnorm):
            pyaws, planes, pconfs = np.append(pyaws,tnorm), np.append(planes,-1), np.append(pconfs,1)
        yaw, yawerr = tyaw, 0
        wRq = geom.RfromYPR(yaw,p,r) # camera orientation (camera to world)
    ### --- THIS IS FOR CHEATING --- ###

    # print pre-homography data to file
    vyaw_err = np.round(np.round(np.mod(tyaw-vyaw,360))) if not np.isnan(vyaw) else np.nan
    vyaw_err = vyaw_err if vyaw_err<180 else 360-vyaw_err
    dbears = np.mod( 180/np.pi*np.arctan2(planeprms[:,0],planeprms[:,2]) , 360 )
    print 'Computed / ground truth cell phone yaw: %.0f / %.0f' % (vyaw,tyaw)
    with open(os.path.join(param['resultsdir'],param['extras_file']),'a') as extras_file:
        print >>extras_file, '\t'.join([timename, '%.0f' % tnorm, '%.0f' % np.round(tyaw), '%.0f' % cyaw, '%.0f' % vyaw, '%.4f'%np.nan, str(vyaw_err)])
        print >>extras_file, '\t'.join([ '%.4f' % 0 for vnorm in vnorms ])
        print >>extras_file, '\t'.join([ '%.0f'   % vnorm  for vnorm  in vnorms  ])
        print >>extras_file, '\t'.join([ '%.0f'   % plane  for plane  in planes  ])
        print >>extras_file, '\t'.join([ '%.0f'   % dplane for dplane in dplanes ])
        print >>extras_file, '\t'.join([ '%.0f'   % dbear  for dbear  in dbears  ])
        print >>extras_file, '\t'.join([ '%.3f' % dnerr  for dnerr  in planeprms[:,4] ])

    # Fill out match information
    nmat = matches['nmat']
    matches['qray'] = geom.normalrows(tp(np.dot(Kqinv,np.append(tp(matches['q2d']),[np.ones(nmat)],0))))
    matches['dray'] = geom.normalrows(tp(np.dot(Kdinv,np.append(tp(matches['d2d']),[np.ones(nmat)],0))))
    matches = match_info(C, Q, matches, dbmatch, wRd)
    matches_start = matches.copy()

    # Solve for query pose using constrained image geometry
    init_matches = initMatches(matches.copy())
    matches['numi'], matches['hconf'] = 0, 0
    runflag, ntry, planechose = param['runflag'], 0, 0
    parameters = ( wRq, wRd, yaw, np.nan, runflag, param['inlier_error'], param['ransac_iter'], 10, True )
    if param['ransac_iter'] == 0:
        matches = init_matches
        matches['numi'], matches['hconf'] == 0, 0
        pose = np.zeros(6)
        pose[3:5] = np.nan
    elif runflag < 10:
        matches, pose = solveGeom(init_matches,parameters,yawerr)
    else:
        ntry = 1
        viter = -np.ones(3)
        parameters = ( wRq, wRd, yaw, np.nan, runflag, param['inlier_error'], param['ransac_iter'], 15, True )
        matches, pose, planechose = solveNorm(C,Q,dbmatch,pyaws,planes,init_matches,parameters,yawerr)
        viter[0] = matches['viter']
        if matches['numi'] == 0 or matches['hconf'] == 0:
            ntry = 2
            parameters = ( wRq, wRd, yaw, np.nan, runflag, 3*param['inlier_error'], param['ransac_iter'], 10, True )
            matches, pose, planechose = solveNorm(C,Q,dbmatch,pyaws,planes,init_matches,parameters,yawerr)
            viter[1] = matches['viter']
        if matches['numi'] == 0 or matches['hconf'] == 0:
            ntry, planechose = 3, 0
            parameters = ( wRq, wRd, yaw, np.nan, 7, 3*param['inlier_error'], param['ransac_iter'], 10, True )
            matches, pose = solveYaw(init_matches,parameters,yawerr)
            viter[2] = matches['viter']
        if matches['numi'] == 0 or matches['hconf'] == 0:
            ntry, planechose = 4, 0

    # save matches to disk
    matches_file = os.path.join(param['resultsdir'],'matches_'+qname+'.pkl')
    matches_out = open(matches_file,'wb')
    pickle.dump(matches,matches_out)
    matches_out.close()
    
    # extract pose parameters
    comp_runflag = matches['runflag']
    tray = pose[:3]
    comp_yaw = pose[3]
    comp_pyaw = pose[4] if runflag>=0 else np.nan
    scale = pose[5] if runflag>=0 else np.nan

    # Get scaled translation for query location
    if np.isnan(scale):
        wRq_pr = geom.YPRfromR(wRq)[1:]
        comp_wRq = geom.RfromYPR(comp_yaw, wRq_pr[0], wRq_pr[1])
        qloc = scalefrom3d(matches, tray, comp_wRq)[0]
    else: # scale calculated in RANSAC loop
        qloc = scale*tray

    # temporarily get yaw error
    qyaw_error = np.round(abs(np.mod(tyaw-comp_yaw,360)))
    qyaw_error = qyaw_error if qyaw_error<180 else abs(qyaw_error-360)

    # compute location errors wrt estimated query locations
    loc_err = ( (qloc[0]-qzx[1])**2 + (qloc[2]-qzx[0])**2 )**0.5
    gps_err = ( (gzx[1] -qzx[1])**2 + (gzx[0] -qzx[0])**2 )**0.5

    # compute the angle difference between T and ground truth translation
    tyaw_error = np.round(abs( 180/np.pi * np.arccos( np.abs(qloc[0]*qzx[1]+qloc[2]*qzx[0]) / (alg.norm([qloc[0],qloc[2]])*alg.norm(qzx)) ) ))

    # compute the plane normal angle error
    nyaw_error = np.nan if np.isnan(comp_pyaw) or np.isnan(tnorm) else np.mod(np.round(abs(comp_pyaw-tnorm)),180)
    nyaw_error = nyaw_error if nyaw_error<90 else abs(nyaw_error-180)

    # write pose estimation results to file
    yaw_err = np.nan
    pose_file = os.path.join(param['resultsdir'],param['pose_file'])
    with open(pose_file,'a') as pf:
        print >>pf, '\t'.join([qname, str(loc_err), str(gps_err), \
            str(tyaw_error), str(qyaw_error), str(nyaw_error), str(matches['numi']), \
            str(matches['numq']), str(matches['nmat']), str(matches['hconf']), \
            str(qloc[0]), str(qloc[2]), str(yaw_err), str(runflag)])

    # print post-homography data to file
    with open(Q.datafile,'a') as df:
        print >>df, ''
        print >>df, '------------------'
        print >>df, ''
        if ntry==1:   print >>df, 'Homography solution using low error threshold with restrictions.'
        elif ntry==2: print >>df, 'Homography solution using high error threshold with restrictions.'
        else:         print >>df, 'Solution not found. Setting T=0.'
        if planechose==0: print >>df, 'Solution formed with unset plane normal.'
        else:                         'Solution chosen with plane normal %d chosen.' % planechose
        print >>df, 'VP yaw | Computed yaw | Actual Yaw | Error : %3.0f | %3.0f | %3.0f | %3.0f' % (vyaw,comp_yaw,tyaw,qyaw_error)
        print >>df, 'Computed Normal | Actual Normal | Error : %3.0f | %3.0f | %3.0f' % (comp_pyaw,tnorm,nyaw_error)
        print >>df, 'Translation   (x|y|z): %.1f | %.1f | %.1f' % (qloc[0],qloc[1],qloc[2])
        print >>df, 'True position (x|-|z): %.1f |  -  | %.1f' % (qzx[1],qzx[0])
        print >>df, 'Angle error | Location error: %.0f | %.1f' % (tyaw_error,loc_err)
        print >>df, 'Number of Inliers | Total matches | Ratio: %d | %d | %.2f' % (matches['numi'],matches['nmat'],np.nan if matches['nmat']==0 else float(matches['numi'])/matches['nmat'])
        print >>df, 'Reprojection error | Homography confidence: %.3f | %.1f' % (matches['rperr'],matches['hconf'])
        print >>df, 'Valid Homographies | Iterations | Ratio: %d | %d | %.3f' % (matches['viter'],matches['niter'],np.nan if matches['niter']==0 else float(matches['viter'])/matches['niter'])
        print >>df, ''
        print >>df, '------------------'
        print >>df, ''
        booleans = [ loc_err<5, loc_err<10, not(5<np.mod(vyaw-tyaw,360)<355), \
                     not(10<np.mod(vyaw-tyaw,360)<350), \
                     not(5<np.mod(comp_yaw-tyaw,360)<355), ntry==1, ntry!=3, \
                     planechose!=0, matches['nmat']!=matches_start['nmat'], \
                     0 if planechose==0 else pconfs[planechose-1]>0, comp_yaw-vyaw ]
        print >>df, '|'.join(['%.0f' % (b) for b in booleans])

    # draw matches
    close = int(loc_err<5) + int(loc_err<10)
    yawclose = int(not(5<np.mod(vyaw-tyaw,360)<355)) + int(not(10<np.mod(vyaw-tyaw,360)<350))
    imgpath = os.path.join( param['resultsdir'] , qname + ';locerr=%.2f' % (loc_err) + ';locPerf_' + str(close) \
        + ';yawPerf_' + str(yawclose) + ';nplanes_' + str(len(pyaws)) + ';plane_' + str(planes[planechose]) + ';try_' + str(ntry) \
        + ';tAng=%.0f' % (tyaw_error) + ';qAng=%.0f' % (qyaw_error) + ';nAng=%.0f' % (nyaw_error) + ';' + dbmatch + '.jpg')
    draw_matches(matches, qimg, dbimg, imgpath)

    print 'Computed yaw / ground truth yaw        : %.0f / %.0f' % (comp_yaw,tyaw)
    if runflag < 10: print 'Computed normal bearing / ground truth : %.0f / %.0f' % (comp_pyaw,tnorm)
    print 'Computed query location relative to db     : %.1f, %.1f, %.1f' % tuple(qloc)
    print 'Ground truth query location relative to db : %.1f,  - , %.1f' % (qzx[1],qzx[0])

    input = (wRq, wRd)

    return qloc, loc_err, matches, input


def initMatches(matches):
    matches['domplane'] = -1
    matches['ifrc'] = 0
    matches['iconf'] = 0
    matches['rperr'] = np.inf
    matches['hconf'] = 0
    matches['numi'] = 0
    matches['viter'] = 0
    matches['niter'] = 0
    matches['runflag'] = -1
    return matches


def solveNorm(C, Q, dbmatch, pyaws, planes, matches, parameters, yawerr):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter, minI, yrestrict = parameters

    if len(pyaws)==0:
        m, p = matches, np.zeros(6)
        m['numi'], m['hconf'], m['viter'], p[3:5] = 0, 0, 0, np.nan
        return m, p, 0
    
    init_nmat = matches['nmat']
    rf = 3 if runflag==10 else 7
    plane_matches, plane_poses, hinls, herrs, hconfs = [], [], np.array([]), np.array([]), np.array([])
    for pyaw,plane in zip(*(pyaws,planes)):
        parameters = ( wRq, wRd, yaw, pyaw, rf, inlerr, rsiter, minI, yrestrict )
        plane_matches.append(match_info(C, Q, matches.copy(), dbmatch, wRd, plane))
        if plane != -1 and plane_matches[-1]['nmat'] < 0.3*init_nmat: # don't allow planes to restrict features beyond this point
            print 'Not using plane because it restricts features too greatly.'
            plane_matches[-1]['nmat'] = 0
            m = plane_matches[-1]
            p = np.zeros(6)
            p[3:5] = np.nan
            plane_poses.append(p)
            hinls = np.append(hinls,m['iconf'])
            herrs = np.append(herrs,m['rperr'])
            hconfs = np.append(hconfs,m['hconf'])
            continue
        print 'Normal yaw set to %.0f degrees' % pyaw
        m, p = solveYaw(plane_matches[-1],parameters,yawerr)
        tray = p[:3]
        wRq_pr = geom.YPRfromR(wRq)[1:]
        comp_wRq = geom.RfromYPR(p[3], wRq_pr[0], wRq_pr[1])
        if np.isnan(p[5]): p[5] = scalefrom3d(m, tray, comp_wRq)[1]
        if abs(p[5])>75:
            m['numi'], m['ifrc'], m['iconf'], m['hconf'], m['rperr'] == 0, 0, 0, 0, np.inf
            p[:] = 0
            p[3:5] = np.nan
        hinls = np.append(hinls,m['iconf'])
        herrs = np.append(herrs,m['rperr'])
        hconfs = np.append(hconfs,m['hconf'])
        plane_matches[-1] = m
        plane_poses.append(p)

    idx = np.argmax(hinls)
    return plane_matches[idx], plane_poses[idx], idx


def solveGeom(matches, parameters, yawerr):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter, minI, yrestrict = parameters
    if runflag < 0: # basic essential matrix
        print 'Basic essential matrix method...'
        runflag = -(1+runflag)
        parameters = ( wRq, wRd, yaw, np.nan, runflag, inlerr, rsiter, yrestrict )
        matches, pose = solveEss(matches,parameters)
    elif runflag < 8: # basic homography
        print 'Basic homography method...'
        if np.isnan(pyaw) and (runflag==1 or runflag==3): runflag -= 1
        parameters = ( wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter, minI, yrestrict )
        matches, pose = solveHom(matches,parameters)
    elif runflag == 8: # varying cell phone yaw outside loop, unknown plane
        print 'Discrete yaw method...'
        parameters = ( wRq, wRd, yaw, pyaw, 3, inlerr, rsiter, minI, yrestrict )
        matches, pose = solveYaw(matches,parameters,yawerr)
    return matches, pose


def solveYaw(matches, parameters, yawerr):
    wRq, wRd, yaw, pyaw, rf, inlerr, rsiter, minI, yrestrict = parameters
    pr = geom.YPRfromR(wRq)[1:]
    dy = 5
    dyawlimit = int(yawerr/dy)*dy
    dyaws = range(-dyawlimit,dyawlimit+dy,dy)
    yaw_matches, yaw_poses, yinls, yerrs, yconfs = [], [], np.array([]), np.array([]), np.array([])
    for dyaw in dyaws:
        yaw_wRq = geom.RfromYPR(yaw+dyaw, pr[0], pr[1])
        runflag = rf - np.isnan(pyaw)
        pyaw_in = pyaw+dyaw if matches['domplane']==-1 else pyaw
        parameters = ( yaw_wRq, wRd, yaw+dyaw, pyaw_in, runflag, inlerr, rsiter, minI, yrestrict )
        yaw_matches.append(matches.copy())
        print 'Query yaw set to %.0f degrees' % (yaw+dyaw)
        m, p = solveHom(yaw_matches[-1],parameters)
        tray = p[:3]
        wRq_pr = geom.YPRfromR(wRq)[1:]
        comp_wRq = geom.RfromYPR(p[3], wRq_pr[0], wRq_pr[1])
        if np.isnan(p[5]): p[5] = scalefrom3d(m, tray, comp_wRq)[1]
        if abs(p[5])>75:
            m['numi'], m['ifrc'], m['iconf'], m['hconf'], m['rperr'] == 0, 0, 0, 0, np.inf
            p[:] = 0
            p[3:5] = np.nan
        yinls = np.append(yinls,m['iconf'])
        yerrs = np.append(yerrs,m['rperr'])
        yconfs = np.append(yconfs,m['hconf'])
        yaw_matches[-1] = m
        yaw_poses.append(p)
    
    ifrc = np.max(yinls)
    mask = yinls > 0.95*ifrc
    yerrs[np.logical_not(mask)] = np.inf
    maskidx = np.argmin(yerrs)
    return yaw_matches[maskidx], yaw_poses[maskidx]


def solveHom(matches, parameters):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter, minI, yrestrict = parameters
    matches, pose = solveHomography.constrainedHomography( matches , \
        wRd , wRq , yaw , pyaw , runflag , maxerr=inlerr , maxiter=rsiter, minI=minI, yrestrict=yrestrict)
    return matches, pose


def solveEss(matches, parameters):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter, yrestrict = parameters
    matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
        wRd , wRq , yaw , runflag , maxerr=inlerr , maxiter=rsiter)
    if matches['numi']==0:
        matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
            wRd , wRq , yaw , runflag , maxerr=3*inlerr , maxiter=rsiter)
    if matches['numi']==0:
        matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
            wRd , wRq , yaw , runflag , maxerr=10*inlerr , maxiter=rsiter)
    return matches, pose


def scalefrom3d(matches, tray, wRq):

    # extract inliers
    imask = matches['imask']
    qray = matches['qray'][imask,:]
    w3d = matches['w3d'][imask,:]
    weights = matches['weight']

    # compute direction of 3d point from query image
    wray = tp(np.dot(wRq,tp(qray)))

    # set up Aq=k equation representing intersection of 2 lines
    numi = np.sum(imask)
    A = [ np.array([[ wray[i,2], -wray[i,0] ], [ tray[2], -tray[0] ]]) for i in range(numi) ]
    k = [ np.array( [ wray[i,2]*w3d[i,0]-wray[i,0]*w3d[i,2] , 0 ] ) for i in range(numi) ]

    # solve for intersection of 2 lines to get query location
    t_int = [ alg.solve(A[i],k[i]) for i in range(numi) ]
    
    # compute the corresponding scale factors attached to y
    idx = int(tray[2]>tray[0])
    scales = [ t_int[i][idx]/tray[2*idx] for i in range(numi) ]

    # find best collection of scale factors
    affprm = [ 2.5 , 0 ] # ransac loop chooses all scale factors within affprm[0]+s*affprm[1] meters of chose scale factor s
    bconf, bnum, bmask = 0, 0, np.bool_(np.zeros(len(scales)))
    for i in range(len(scales)):
        s = scales[i]
        mask = np.abs(s-scales) < affprm[0] + affprm[1] * np.abs(s)
        if np.sum(weights[mask]) > bconf : bconf, bnum, bmask = np.sum(weights[mask]), np.sum(mask), mask
    tdist = np.mean(np.compress(bmask,scales))
    t = tdist*tray

    print '%.0f / %.0f matches used to compute scale.' % (bnum,numi)
    
    return t, tdist
