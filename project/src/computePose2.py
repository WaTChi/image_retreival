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


def highresSift(C, Q, dbmatch):

    # timing info
    start = time.time()

    # set sift paths
    qname = Q.name
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
    matches['ddep'] = np.zeros(nmat)
    matches['plane'] = np.int_( -1 * np.ones(nmat) )
    matches['weight'] = np.zeros(nmat)
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


def find_domplanes(C, Q, dbmatch, Kdinv, wRd):
    print 'Finding dominant planes...'
    if C.QUERY != 'oakland1': return ( np.zeros(0), np.zeros(0), np.zeros((0,5)) )
    planedata = np.asarray( Image.open( os.path.join(C.hiresdir,dbmatch+'-planes.png') ) )
    h,w = planedata.shape
    maxplane = np.max(planedata[planedata!=255])
    plane_count = np.array( [ float(np.sum(np.reshape(planedata==k+1,-1)))/(h*w) for k in range(maxplane) ] )
    planes = 1 + np.nonzero(plane_count>0.1)[0]
    sizes = plane_count[planes-1]
    if len(planes) == 0: return planes, sizes, np.zeros((0,5),np.float)
    print 'Computing dominant plane normals...'
    start = time.time()
    params = np.array( [ planefrom3d(C, Q, dbmatch, planes[i], Kdinv, wRd) for i in range(len(planes)) ] )
    print 'Computing plane normals took %.1f seconds.' % (time.time()-start)
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
    p3d_norm = geom.normalrows(p3d)

    # extract least squared normal vector
    neig = alg.eig(np.dot(tp(p3d_norm),p3d_norm))
    n = neig[1][:3,np.argmin(neig[0])]

    # restrict normal vector to be perpendicular to gravity
    g = np.array([0,1,0])
    n = geom.normalrows( n - np.inner(n,g) * g )
    
    # get plane distance and flip n if needed
    pd = -np.mean(np.inner(p3d[:,:3],n))
    if pd<0: n, pd = -n, -pd

    # get average point error in meters
    err = np.mean(np.abs(np.inner(p3d,np.append(n,pd))))

    return np.append(n,[pd,err])


def combine_planes(runflag,vnorms,vpconfs,dplanes,dsizes,dprms):
    pyaws, planes, pconfs = np.zeros(0), np.zeros(0,np.int), np.zeros(0)
    rmidx = np.nonzero(vpconfs<0.51)[0]
    vnorms, vpconfs = np.delete(vnorms,rmidx), np.delete(vpconfs,rmidx)
    while len(vnorms) > 0 and runflag > 6:
        vyaw = vnorms[0]
        dyaws = 180/np.pi*np.arctan2(dprms[:,0],dprms[:,2])
        diffs = np.mod(vyaw-dyaws,360)
        diffs[diffs>180] = 360-diffs[diffs>180]
        mask = diffs < 20
        plane, psize = dplanes[mask], dsizes[mask]
        pyaws, pconfs = np.append(pyaws,vyaw), np.append(pconfs,vpconfs[0])
        if len(plane)==0: planes = np.append(planes,-1)
        else: planes = np.append(planes,plane[np.argmax(psize)])
        vnorms, vpconfs = np.delete(vnorms,0), np.delete(vpconfs,0)
        rmidx = np.nonzero(mask)[0]
        if len(rmidx)>0: dplanes, dsizes, dprms = np.delete(dplanes,rmidx), np.delete(dsizes,rmidx), np.delete(dprms,rmidx,0)
    while len(dplanes) > 0:
        dyaw = np.mod(180/np.pi*np.arctan2(dprms[0,0],dprms[0,2]),360)
        if dprms[0,-1] < 5:
            pyaws, planes, pconfs = np.append(pyaws,dyaw), np.append(planes,dplanes[0]), np.append(pconfs,0)
        dplanes, dsizes, dprms = np.delete(dplanes,0), np.delete(dsizes,0), np.delete(dprms,0,0)
    return pyaws, planes, pconfs


def estimate_pose(C, Q, dbmatch, gtStatus=None):

    # get pose parameters
    param = C.pose_param
    runflag = param['runflag']
    print param

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

    # get high res sift rematch
    matches = highresSift(C, Q, dbmatch)

    # Get estimated ground truth query location and normal direction
    tlat, tlon, tnorm, tyaw = getGTpose(C, Q)
    qzx = geom.lltom(olat,olon,tlat,tlon)

    # get query yaw and plane yaw from vanishing point anaylsis
    vyaw, vyawconf, vps, vpcenters, vnorms, vpconfs, nqvps = vp_analysis.getQNyaws(C, Q, qimg, dbimg, qsource)

    # get dominant planes
    dplanes, psizes, planeprms = find_domplanes(C, Q, dbmatch, Kdinv, wRd)
    print 'Dominant planes: ' + str(dplanes)
    
    # match vanishing point planes to database planes
    pyaws, planes, pconfs = combine_planes(runflag,vnorms,vpconfs,dplanes,psizes,planeprms)

    # Set Kq, wRq
    wx,wy = qsource.image.size
    fov = qsource.fov
    Kq = geom.cameramat(wx, wy, fov)
    Kqinv = alg.inv(Kq)
    cyaw,p,r = qsource.yaw, qsource.pitch, qsource.roll # cyaw - cell phone yaw

    # Set yaw value to be used
    if runflag > 6: # vanishing point methods
        if vyawconf < 0.45: runflag, yaw, yawerr = 6, cyaw, 30
        else: runflag, yaw, yawerr = 5, vyaw, 0
    else: yaw, yawerr = cyaw, 30 # set cell phone yaw to use, plane normal
#    runflag, yaw, yawerr = 5, tyaw, 0
    wRq = geom.RfromYPR(yaw,p,r) # camera orientation (camera to world)

    ### --- THIS IS FOR CHEATING --- ###
    if param['cheat']:
        if np.isnan(tnorm):
            pyaws, planes, pconfs = np.zeros(0), np.zeros(0,np.int), np.zeros(0)
        else:
            cheatdiff = np.mod(pyaws-tnorm,360)
            cheatdiff[cheatdiff>180] = 360-cheatdiff[cheatdiff>180]
            cheatmask = cheatdiff<20
            cheatyaws, cheatplanes,cheatconfs = pyaws[cheatmask], planes[cheatmask], pconfs[cheatmask]
            if len(cheatyaws) == 0:
                pyaws, planes, pconfs = np.array([tnorm]), -1*np.ones(1,np.int), np.zeros(1)
            else:
                cheatidx = np.argmax(cheatconfs)
                pyaws, planes, pconfs = np.array([tnorm]), np.array([cheatplanes[cheatidx]]), np.array([cheatconfs[cheatidx]])
        runflag, yaw, yawerr = 5, tyaw, 0
        wRq = geom.RfromYPR(yaw,p,r) # camera orientation (camera to world)
    ### --- THIS IS FOR CHEATING --- ###

    # print pre-homography data to file
    vyaw_err = int(np.round(np.mod(tyaw-vyaw,360))) if not np.isnan(vyaw) else np.nan
    vyaw_err = vyaw_err if vyaw_err<180 else 360-vyaw_err
    dbears = np.mod( 180/np.pi*np.arctan2(planeprms[:,0],planeprms[:,2]) , 360 )
    print 'Computed / ground truth cell phone yaw: %.0f / %.0f' % (vyaw,tyaw)
    with open(param['extras_file'],'a') as extras_file:
        print >>extras_file, '\t'.join([timename, '%.0f' % tnorm, str(int(tyaw)), '%.0f' % cyaw, '%.0f' % vyaw, '%.4f'%vyawconf, str(vyaw_err)])
        print >>extras_file, '\t'.join([ '%.4f' % vpconf for vpconf in vpconfs ])
        print >>extras_file, '\t'.join([ '%d'   % vnorm  for vnorm  in vnorms  ])
        print >>extras_file, '\t'.join([ '%d'   % plane  for plane  in planes  ])
        print >>extras_file, '\t'.join([ '%d'   % dplane for dplane in dplanes ])
        print >>extras_file, '\t'.join([ '%d'   % dbear  for dbear  in dbears  ])
        print >>extras_file, '\t'.join([ '%.3f' % dnerr  for dnerr  in planeprms[:,4] ])

    # Fill out match information
    nmat = matches['nmat']
    matches['qray'] = geom.normalrows(tp(np.dot(Kqinv,np.append(tp(matches['q2d']),[np.ones(nmat)],0))))
    matches['dray'] = geom.normalrows(tp(np.dot(Kdinv,np.append(tp(matches['d2d']),[np.ones(nmat)],0))))
    matches = match_info(C, Q, matches, dbmatch, wRd)

    # Solve for query pose using constrained image geometry
    runflag = param['runflag']
    parameters = ( wRq, wRd, yaw, np.nan, runflag, param['inlier_error'], param['ransac_iter'] )
    if runflag < 5: matches, pose = solveGeom(matches,parameters,yawerr)
    else: matches, pose = solveNorm(C,Q,dbmatch,pyaws,planes,pconfs,matches,parameters,yawerr,qzx,tyaw)
    
    # extract pose parameters
    tray = pose[:3]
    comp_yaw = pose[3]
    comp_pyaw = pose[4] if runflag<10 else np.nan

    # Get scaled translation for query location
    wRq_pr = geom.YPRfromR(wRq)[1:]
    comp_wRq = geom.RfromYPR(comp_yaw, wRq_pr[0], wRq_pr[1])
    qloc = scalefrom3d(matches, tray, comp_wRq)[0]

    # temporarily get yaw error
    qyaw_error = int(abs(np.mod(tyaw-comp_yaw,360)))
    qyaw_error = qyaw_error if qyaw_error<180 else abs(qyaw_error-360)

    # compute location errors wrt estimated query locations
    loc_err = ( (qloc[0]-qzx[1])**2 + (qloc[2]-qzx[0])**2 )**0.5
    gps_err = ( (gzx[1] -qzx[1])**2 + (gzx[0] -qzx[0])**2 )**0.5

    # compute the angle difference between T and ground truth translation
    tyaw_error = int(abs( 180/np.pi * np.arccos( np.abs(qloc[0]*qzx[1]+qloc[2]*qzx[0]) / (alg.norm([qloc[0],qloc[2]])*alg.norm(qzx)) ) ))

    # compute the plane normal angle error
    nyaw_error = np.nan if np.isnan(comp_pyaw) or np.isnan(tnorm) else np.mod(int(abs(comp_pyaw-tnorm)),180)
    nyaw_error = nyaw_error if nyaw_error<90 else abs(nyaw_error-180)

    # write pose estimation results to file
    yaw_err = np.nan
    with open(param['pose_file'],'a') as pose_file:
        print >>pose_file, '\t'.join([qname, str(loc_err), str(gps_err), \
            str(tyaw_error), str(qyaw_error), str(nyaw_error), str(matches['numi']), \
            str(matches['numq']), str(matches['nmat']), str(matches['hconf']), \
            str(qloc[0]), str(qloc[2]), str(yaw_err), str(runflag)])

    # draw matches
    close = loc_err < 10
    imgpath = os.path.join( param['resultsdir'] , qname + ';nplanes=' + str(nqvps) + ';gt' + str(gtStatus) + ';' + 'close' + str(close) + ';locerr=' + str(loc_err) + ';tAng=' + str(tyaw_error) + ';qAng=' + str(qyaw_error) + ';nAng=' + str(nyaw_error) + ';' + dbmatch + '.jpg')
    draw_matches(matches, qimg, dbimg, imgpath)

    print 'Computed yaw / ground truth yaw        : %.0f / %.0f' % (comp_yaw,tyaw)
    if runflag < 10: print 'Computed normal bearing / ground truth : %.0f / %.0f' % (comp_pyaw,tnorm)
    print 'Computed query location relative to db     : %.1f, %.1f, %.1f' % tuple(qloc)
    print 'Ground truth query location relative to db : %.1f,  - , %.1f' % (qzx[1],qzx[0])

    input = (wRq, wRd)

    return qloc, loc_err, matches, input


def solveNorm(C, Q, dbmatch, pyaws, planes, pconfs, matches, parameters, yawerr, qzx, tyaw):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter = parameters

    ### Print planes and location of query ###
    planeprint = ['%d'%tyaw]
    for py in pyaws: planeprint.append('%d'%py)
    with open(Q.datafile,'a') as df:
        print >>df, '\t'.join(planeprint)
        print >>df, '\t'.join(['%d'%qzx[1],'-','%d'%qzx[0]])
        print >>df, ''
    ### Print planes and location of query ###

#    planemask = planeprms[:,4] < 3
#    planes = planes[planemask]
#    pyaws = np.mod( 180 / np.pi * np.arctan2(planeprms[planemask,0],planeprms[planemask,2]) , 360 )
    plane_matches, plane_poses, hconfs = [], [], []
    for pyaw,plane,pconf in zip(*(pyaws,planes,pconfs)):
        parameters = ( wRq, wRd, yaw, pyaw, 4, inlerr, rsiter )
        plane_matches.append(match_info(C, Q, matches.copy(), dbmatch, wRd, plane))
        print 'Normal yaw set to %d degrees' % pyaw
        m, p = solveYaw(plane_matches[-1],parameters,yawerr)
#        if runflag == 5: m, p = solveHom(plane_matches[-1],parameters) # cell yaw set
#        else: m, p = solveYaw(plane_matches[-1],parameters) # runflag == 6, solve cell yaw
        plane_matches[-1] = m
        plane_poses.append(p)
        weight = 1.6 if pconf>0 and plane!=-1 else 1.3
        hconfs.append(weight*m['hconf']) # if m['numi']>=10 else hconfs.append(0)

        ### Get scale, error, and print results ###
        tray = p[:3]
        wRq_pr = geom.YPRfromR(wRq)[1:]
        comp_wRq = geom.RfromYPR(p[3], wRq_pr[0], wRq_pr[1])
        qloc = scalefrom3d(m, tray, comp_wRq)[0]
        qyaw_error = int(abs(np.mod(tyaw-p[3],360)))
        qyaw_error = qyaw_error if qyaw_error<180 else abs(qyaw_error-360)
        tyaw_error = int(abs( 180/np.pi * np.arccos( np.abs(qloc[0]*qzx[1]+qloc[2]*qzx[0]) / (alg.norm([qloc[0],qloc[2]])*alg.norm(qzx)) ) ))
        rel_norm = np.mod(int(p[3]+180-p[4]),360)
        rel_norm = rel_norm if rel_norm < 180 else rel_norm-360
        with open(Q.datafile,'a') as df:
            print >>df, '\t'.join(['%d'%rel_norm,'%d'%qyaw_error,'%d'%tyaw_error])
            print >>df, '\t'.join(['%d'%qloc[0],'%d'%qloc[1],'%d'%qloc[2]])
            print >>df, '\t'.join(['%.1f'%abs(qloc[0]-qzx[1]),'%.1f'%abs(2-qloc[1]),'%.1f'%abs(qloc[2]-qzx[0])])
            print >>df, '\t'.join(['%.1f'%m['hconf'],'%.5f'%m['rperr']])
            print >>df, '\t'.join(['%d'%(100.*m['numi']/m['numq']),'%d'%m['numi'],'%d'%m['numq']])
            print >>df, ''
        ### Get scale, error, and print results ###


    parameters = ( wRq, wRd, yaw, np.nan, 4, inlerr, rsiter )
    plane_matches.append(matches.copy())
    print 'Normal yaw unset'
    m, p = solveYaw(plane_matches[-1],parameters,yawerr)
#    if runflag == 5: m, p = solveHom(plane_matches[-1],parameters) # cell yaw set
#    else: m, p = solveYaw(plane_matches[-1],parameters) # runflag == 6, solve cell yaw
    plane_matches[-1] = m
    plane_poses.append(p)
    hconfs.append(m['hconf'])

    ### Get scale, error, and print results ###
    tray = p[:3]
    wRq_pr = geom.YPRfromR(wRq)[1:]
    comp_wRq = geom.RfromYPR(p[3], wRq_pr[0], wRq_pr[1])
    qloc = scalefrom3d(m, tray, comp_wRq)[0]
    qyaw_error = int(abs(np.mod(tyaw-p[3],360)))
    qyaw_error = qyaw_error if qyaw_error<180 else abs(qyaw_error-360)
    tyaw_error = int(abs( 180/np.pi * np.arccos( np.abs(qloc[0]*qzx[1]+qloc[2]*qzx[0]) / (alg.norm([qloc[0],qloc[2]])*alg.norm(qzx)) ) ))
    rel_norm = np.mod(int(p[4]+180-p[3]),360)
    rel_norm = rel_norm if rel_norm < 180 else rel_norm-360
    with open(Q.datafile,'a') as df:
        print >>df, '\t'.join(['-','-','-'])
        print >>df, ''
        print >>df, '\t'.join(['%d'%rel_norm,'%d'%qyaw_error,'%d'%tyaw_error])
        print >>df, '\t'.join(['%d'%qloc[0],'%d'%qloc[1],'%d'%qloc[2]])
        print >>df, '\t'.join(['%.1f'%abs(qloc[0]-qzx[1]),'%.1f'%abs(2-qloc[1]),'%.1f'%abs(qloc[2]-qzx[0])])
        print >>df, '\t'.join(['%.1f'%m['hconf'],'%.5f'%m['rperr']])
        print >>df, '\t'.join(['%d'%(100.*m['numi']/m['numq']),'%d'%m['numi'],'%d'%m['numq']])
        print >>df, ''
    ### Get scale, error, and print results ###

    bestidx = np.argmax(hconfs)
    return plane_matches[bestidx], plane_poses[bestidx]


def solveGeom(matches, parameters, yawerr):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter = parameters
    if runflag < 0: # basic essential matrix
        print 'Basic essential matrix method...'
        runflag = -(1+runflag)
        parameters = ( wRq, wRd, yaw, np.nan, runflag, inlerr, rsiter )
        matches, pose = solveEss(matches,parameters)
    elif runflag < 4: # basic homography
        print 'Basic homography method...'
        if np.isnan(pyaw) and (runflag==1 or runflag==3): runflag -= 1
        parameters = ( wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter )
        matches, pose = solveHom(matches,parameters)
    elif runflag == 4: # varying cell phone yaw outside loop, unknown plane
        print 'Discrete yaw method...'
        parameters = ( wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter )
        matches, pose = solveYaw(matches,parameters,yawerr)
    return matches, pose


def solveYaw(matches, parameters, yawerr):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter = parameters
    pr = geom.YPRfromR(wRq)[1:]
    dy = 3
    dyawlimit = int(yawerr/dy)*dy
    dyaws = range(-dyawlimit,dyawlimit+dy,dy)
    rsiter = int( rsiter / len(dyaws)**0.5 )
    yaw_matches, yaw_poses, yconfs = [], [], []
    for dyaw in dyaws:
        yaw_wRq = geom.RfromYPR(yaw+dyaw, pr[0], pr[1])
        runflag = 3 - np.isnan(pyaw)
        parameters = ( yaw_wRq, wRd, yaw+dyaw, pyaw, runflag, inlerr, rsiter )
        yaw_matches.append(matches.copy())
        print 'Query yaw set to %d degrees' % (yaw+dyaw)
        m, p = solveHom(yaw_matches[-1],parameters)
        yaw_matches[-1] = m
        yaw_poses.append(p)
        yconfs.append(m['hconf'])
    bestidx = np.argmax(yconfs)
    return yaw_matches[bestidx], yaw_poses[bestidx]


def solveHom(matches, parameters):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter = parameters
    matches, pose = solveHomography.constrainedHomography( matches , \
        wRd , wRq , yaw , pyaw , runflag , maxerr=inlerr , maxiter=rsiter)
    if matches['numi']==0 or geom.vecnorm(pose[:3])>50:
        matches, pose = solveHomography.constrainedHomography( matches , \
            wRd , wRq , yaw , pyaw , runflag , maxerr=2*inlerr , maxiter=rsiter)
    if matches['numi']==0 or geom.vecnorm(pose[:3])>50:
        matches, pose = solveHomography.constrainedHomography( matches , \
            wRd , wRq , yaw , pyaw , runflag , maxerr=5*inlerr , maxiter=rsiter)
    return matches, pose


def solveEss(matches, parameters):
    wRq, wRd, yaw, pyaw, runflag, inlerr, rsiter = parameters
    matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
        wRd , wRq , yaw , runflag , maxerr=inlerr , maxiter=rsiter)
    if matches['numi']==0 or geom.vecnorm(pose[:3])>50:
        matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
            wRd , wRq , yaw , runflag , maxerr=2*inlerr , maxiter=rsiter)
    if matches['numi']==0 or geom.vecnorm(pose[:3])>50:
        matches, pose = solveEssMatrix.constrainedEssMatrix( matches , \
            wRd , wRq , yaw , runflag , maxerr=5*inlerr , maxiter=rsiter)
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

    print '%d / %d matches used to compute scale.' % (bnum,numi)
    
    return t, tdist