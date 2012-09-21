import geom
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="aaronh"
__date__ ="Nov 29, 2011"

import time
import numpy as np
from numpy import transpose as tp
import numpy.linalg as alg
import numpy.random as rnd
import scipy.optimize as opt
import render_tags
import os


#####                                                                       #####
#####  Main function for getting qyaw, nyaw from vanishing points analysis  #####
#####                                                                       #####
def getQNyaws(C, Q, qimg, dimg, qsource, tyaw):

    print 'Retrieving cell phone yaw and normal bearings from vanishing points...'
    start = time.time()

    # all VPs within the following degree threshold are considered the same VP
    vp_threshold = 5 # should be similar to the yaw error desired

    # get building face horizontal vp from database image(s) ; compute nyaw
    vps, vpcenters, vnorms, vpconfs, ndvps = VPNfromDatabase(C, Q, dimg, vp_threshold)

    # match vp from query to vp above to compute qyaw
    if len(vpconfs)==0: vyaw, yawconf = np.nan, 0
    else: vyaw, yawconf, vps, vnorms, vpconfs, nqvps = VPQfromQuery(C, Q, qimg, qsource, vps, vnorms, vpconfs, vp_threshold, tyaw)

    print 'Vanishing point analysis took %.1f seconds.' % (time.time()-start)
    if len(vpconfs)==0: print 'Database vanishing point analysis failed.'
    elif np.isnan(vyaw): print 'Query vanishing point analyis failed.'
    else:
        print 'Computed cell phone yaw / confidence  : %d / %.2f' % (vyaw,yawconf)
        print 'Computed normal bearings: ' + str([int(round(vn)) for vn in vnorms])


    return vyaw, vnorms


def VPNfromDatabase(C, Q, dimg, vp_threshold):

    main_bias, off_bias = 1, 0

    if off_bias == 0:

        dname = os.path.basename(dimg)
        himg, dinfo, dpath = os.path.join(C.hiresdir, dname[:-4] + '.jpg'), \
                             os.path.join(C.hiresdir, dname[:-4] + '.info'), \
                             os.path.join(C.hiresdir, dname[:-4] + '.lsd')
        dsource = render_tags.EarthmineImageInfo(himg,dinfo)
        Kd, wRd = viewparam(dsource,np.nan)
        dmid, deqs, dlen = LfromLSD(dpath,himg,Kd)
        dvps, dcon, dcent, dseeds = VPfromSeeds(dmid, deqs, dlen, wRd, vp_threshold)
        vps, conf, cent = tp(np.dot(wRd,tp(dvps))), dcon, dcent
        nvps = len(conf)
        if nvps == 0: return np.zeros((0,3)), np.zeros((0,3)), np.zeros(0), np.zeros(0) # return if no vanishing points
    
    else:

        # get 3 database images
        dname = os.path.basename(dimg)
        view = int(dname[-6:-4])
        if view < 6: # right side of street
            limg, linfo, lpath = os.path.join(C.hiresdir, dname[:-6] + '02.jpg'), \
                                 os.path.join(C.hiresdir, dname[:-6] + '02.info'), \
                                 os.path.join(C.hiresdir, 'lsd', dname[:-6] + '02.lsd')
            cimg, cinfo, cpath = os.path.join(C.hiresdir, dname[:-6] + '03.jpg'), \
                                 os.path.join(C.hiresdir, dname[:-6] + '03.info'), \
                                 os.path.join(C.hiresdir, 'lsd', dname[:-6] + '03.lsd')
            rimg, rinfo, rpath = os.path.join(C.hiresdir, dname[:-6] + '04.jpg'), \
                                 os.path.join(C.hiresdir, dname[:-6] + '04.info'), \
                                 os.path.join(C.hiresdir, 'lsd', dname[:-6] + '04.lsd')
        else: # left side of street
            limg, linfo, lpath = os.path.join(C.hiresdir, dname[:-6] + '08.jpg'), \
                                 os.path.join(C.hiresdir, dname[:-6] + '08.info'), \
                                 os.path.join(C.hiresdir, 'lsd', dname[:-6] + '08.lsd')
            cimg, cinfo, cpath = os.path.join(C.hiresdir, dname[:-6] + '09.jpg'), \
                                 os.path.join(C.hiresdir, dname[:-6] + '09.info'), \
                                 os.path.join(C.hiresdir, 'lsd', dname[:-6] + '09.lsd')
            rimg, rinfo, rpath = os.path.join(C.hiresdir, dname[:-6] + '10.jpg'), \
                                 os.path.join(C.hiresdir, dname[:-6] + '10.info'), \
                                 os.path.join(C.hiresdir, 'lsd', dname[:-6] + '10.lsd')
            lsource = render_tags.EarthmineImageInfo(limg, linfo)
            csource = render_tags.EarthmineImageInfo(cimg, cinfo)
            rsource = render_tags.EarthmineImageInfo(rimg, rinfo)

        # extract view parameters
        Kl, wRl = viewparam(lsource,np.nan)
        Kc, wRc = viewparam(csource,np.nan)
        Kr, wRr = viewparam(rsource,np.nan)
        
        # get lines for each database image; image frame equations and segment lengths
        lmid, leqs, llen = LfromLSD(lpath, limg, Kl)
        cmid, ceqs, clen = LfromLSD(cpath, cimg, Kc)
        rmid, reqs, rlen = LfromLSD(rpath, rimg, Kr)

        # get candidate vanishing points from lines
        lvps, lcon, lcent, lseeds = VPfromSeeds(lmid, leqs, llen, wRl, vp_threshold)
        cvps, ccon, ccent, cseeds = VPfromSeeds(cmid, ceqs, clen, wRc, vp_threshold)
        rvps, rcon, rcent, rseeds = VPfromSeeds(rmid, reqs, rlen, wRr, vp_threshold)

        #####  combine candidate vanishing points and into an estimate of   #####
        #####  the building faces' horizontal vanishing points and normals  #####

        # increase the confidence of vps from the matched view and
        if    view==2 or view==8  : lcon, ccon, rcon, ccent, rcent, ndvps, seedlens = main_bias*lcon, off_bias*ccon, off_bias*rcon, 0*ccent, 0*rcent, len(lvps), lseeds
        elif  view==3 or view==9  : lcon, ccon, rcon, lcent, rcent, ndvps, seedlens = off_bias*lcon, main_bias*ccon, off_bias*rcon, 0*lcent, 0*rcent, len(cvps), cseeds
        elif  view==4 or view==10 : lcon, ccon, rcon, lcent, ccent, ndvps, seedlens = off_bias*lcon, off_bias*ccon, main_bias*rcon, 0*lcent, 0*ccent, len(rvps), rseeds

        # map the vanishing points to the world frame (EDN - east/down/north) and combine all vps
        lvps, cvps, rvps = tp(np.dot(wRl,tp(lvps))), tp(np.dot(wRc,tp(cvps))), tp(np.dot(wRr,tp(rvps)))
        vps, conf, cent = np.concatenate( (lvps,cvps,rvps) , 0 ), np.concatenate((lcon,ccon,rcon)), np.concatenate( (lcent,ccent,rcent) , 0 )
        nvps = len(conf)
        if nvps == 0: return np.zeros((0,3)), np.zeros((0,3)), np.zeros(0), np.zeros(0) # return if no vanishing points

    # get normals and remove vanishing points indicating more than a ~18 degree incline
    normals = np.cross(vps,[0,1,0])
    mask = geom.vecnorm(normals) > 0.95
    vps, cent, normals, conf = vps[mask,:], cent[mask,:], geom.normalrows(normals[mask,:]), conf[mask]
    nvps = len(conf)

    # sort vanishing points by confidence
    sort = np.argsort(conf)
    vps, cent, conf = vps[sort[::-1],:], cent[sort[::-1],:], conf[sort[::-1]]

    # combine all vanishing points
    minconf = 0.2  # average 20% of line length in each image OR 50% of line length in retrieved image
    bvps, bcenters, bnorms, bconfs = np.zeros((0,3)), np.zeros((0,3)), np.zeros(0), np.zeros(0)
    while len(conf)!=0:
        vp = vps[0,:]
        mask = np.inner(vps,vp) > np.cos(vp_threshold*np.pi/180)
        c = np.sum(conf[mask])/(2*off_bias+main_bias)
        if c > minconf:
            vp = geom.largestSingVector(geom.vecmul(vps[mask,:],conf[mask]))
            if np.inner(vps[0,:],vp) < 0: vp = -vp
            normal = np.cross(vp,[0,1,0])
            nyaw = np.mod( 180/np.pi * np.arctan2(normal[0],normal[2]) , 360 )
            bvps = np.concatenate( (bvps,[vp]) , 0 )
            bnorms, bconfs = np.append(bnorms,nyaw), np.append(bconfs,c)
            centmask = np.logical_and(mask,cent[:,2]!=0)
            center = np.mean(cent[centmask,:],0)
            bcenters = np.concatenate( (bcenters,[center]) , 0 )
            keep = np.logical_not(mask)
            vps, conf, cent = vps[keep,:], conf[keep], cent[keep,:]
        else:
            vps, conf, cent = np.delete(vps,0,0), np.delete(conf,0), np.delete(cent,0,0)

    # sort best vanishing points by confidence
    if len(bconfs) == 0: return bvps, bcenters, bnorms, bconfs
    sort = np.argsort(bconfs)
    bvps, bcenters, bnorms, bconfs = bvps[sort[::-1],:], bcenters[sort[::-1],:], bnorms[sort[::-1]], bconfs[sort[::-1]]

    return bvps, bcenters, bnorms, bconfs, nvps


def VPQfromQuery(C, Q, qimg, qsource, vps, vnorms, vpconfs, vp_threshold, tyaw):

    # get query vanishing points
    qname = os.path.basename(qimg)
    qpath = os.path.join(C.querydir, 'hires', 'lsd', qname[:-4] + '.lsd')
    Kq, wRq = viewparam(qsource,tyaw)
    qmid, qleq, qlen = LfromLSD(qpath, qimg, Kq)
    qvps, conf, qcent, seedlens = VPfromSeeds(qmid, qleq, qlen, wRq, vp_threshold)
    nqvps = len(conf)

    #####  combine candidate vanishing points and vp from db   #####
    #####  into an estimate of the true query yaw orientation  #####

    # map vanishing points to world frame
    qvps = tp(np.dot(wRq,tp(qvps)))
    
    # align vanishing points based on normal and compute normals
    qnorms = geom.normalrows(np.cross(qvps,[0,1,0]))
    for i in range(len(conf)):
        if np.dot(tp(wRq),qnorms[i,:])[2] > 0:
            qnorms[i,:] *= -1
            qvps[i,:] *= -1

    # find optimal alignment of vanishing points
    cyaw = geom.YPRfromR(wRq)[0] # cell phone yaw
    byaw, bconf, bvps, bnorms, bvpconfs, nvps = np.nan, 0, vps, vnorms, vpconfs, len(vpconfs)

#    print '------------------------'
#    print vpconfs
#    print np.mod( vnorms , 360)
#    print conf
#    print np.mod( 180/np.pi * np.arctan2(qnorms[:,0],qnorms[:,2]) , 360 )
#    print '------------------------'
    qnormyaws = 180/np.pi * np.arctan2(qnorms[:,0],qnorms[:,2])
    for i in range(len(vpconfs)):
        for j in range(len(conf)):
            # compute relative yaw change
            vnormyaw = vnorms[i] #180/np.pi * np.arctan2(vnorms[i,0],vnorms[i,2])
            qnormyaw = qnormyaws[j]
            dyaw = vnormyaw - qnormyaw
            dyaw = dyaw if dyaw<180 else dyaw-360
            if abs(dyaw) > 50: continue # skip if the yaw change is too great
            # apply relative yaw change
            dR = geom.RfromYPR(dyaw,0,0)
            dvps, dnorms = tp(np.dot(dR,tp(qvps))), tp(np.dot(dR,tp(qnorms)))
            # get list of matching vanishing points
            dbidx, qidx, weights = np.zeros(0,np.int), np.zeros(0,np.int), np.zeros(0)
            # Gather lise of aligned vanishing points
            for k in range(len(vpconfs)):
                vpalign = np.inner(dvps,vps[k,:])
                alignidx = np.argmax(vpalign)
                if vpalign[alignidx] < np.cos(np.pi/180*2*vp_threshold): continue
                dbidx, qidx = np.append(dbidx,k), np.append(qidx,alignidx)
                weights = np.append(weights,conf[alignidx]*vpconfs[k])
            # Optimize for the yaw change
            yawconf = np.sum(weights)
            if yawconf <= bconf: continue
            dyaws = np.mod(vnorms[dbidx]-qnormyaws[qidx],360)
            if dyaws[0] < 90: dyaws[dyaws>270] = dyaws[dyaws>270]-360
            elif dyaws[0] > 270: dyaws[dyaws<90] = dyaws[dyaws<90]+360
            dyaw = np.sum(weights*dyaws) / yawconf
            byaw, bconf, bvpconfs, nvps = np.mod(cyaw+dyaw,360), yawconf, np.ones(len(weights)), len(weights)
            bnorms = np.mod( qnormyaws[qidx] + dyaw , 360 )

    return byaw, bconf, bvps, bnorms, bvpconfs, nqvps

def alignVPcost(dyaw,vpmat1,vpmat2,weights):
    dR = geom.RfromYPR(dyaw[0],0,0)
    matchdot = np.sum( vpmat1*tp(np.dot(dR,tp(vpmat2))) , 1 )
    yaw_weight = 1 # np.cos(np.pi/180*dyaw[0])
    return -yaw_weight*np.sum(weights*matchdot)

def viewparam(source,tyaw):
    # camera calibration matrix
    Kcal = geom.cameramat(source.image.size[0], source.image.size[1], source.fov)
    # camera orientation (camera to world)
    yaw = source.yaw if np.isnan(tyaw) else tyaw
    Rot = geom.RfromYPR(yaw, source.pitch, source.roll)
    return Kcal, Rot


def VPfromSeeds(midpts, lineqs, lengths, Rot, tol):

    # Generate seeds from rotation matrix
    total_len = np.sum(lengths)
    seed_tol = 2*tol
    yaw, Rpr = geom.YfromR(Rot)
    angles = np.arange(0,180,3)
    nvps = len(angles)
    vps, vlens = np.zeros((nvps,3)), np.zeros(nvps)
    for i in xrange(nvps):
        vps[i,:] = np.dot( tp(Rpr) , [np.sin(np.pi/180*angles[i]),0,np.cos(np.pi/180*angles[i])] )

    # Iterate through each line and assign its weight to top N seeds
    nseeds = 3
    for i in xrange(len(lengths)):
        midpt, lineq, length = midpts[i,:], lineqs[i,:], lengths[i]
        line_bear = np.arctan(lineq[0]/lineq[1])
        vp_bear = np.arctan( (midpt[1]-vps[:,1]/vps[:,2]) / (vps[:,0]/vps[:,2]-midpt[0]) )
        dbear = np.mod(line_bear-vp_bear,np.pi)
        dbear = dbear + (dbear>np.pi/2) * (np.pi-2*dbear)
        vpdist = geom.vecnorm(geom.vecsub(geom.vecdiv(vps,vps[:,2],0),midpt,1))
        mask = vpdist < length/2
        dbear[mask] = np.pi
        minidx = np.argsort(dbear)[:nseeds]
        vlens[minidx] += length/total_len

    # Pick true vanishing points from seeds
#    if True:
#        file = '/media/DATAPART2/ah/pose_runs/tmp.txt'
#        open(file,'w').close()
#        with open(file,'a') as f:
#            for tmp in vlens: print >>f, '%.3f' % tmp
    seedlens = vlens
    neighbors = np.amax([np.roll(vlens,-3),np.roll(vlens,-2),np.roll(vlens,-1), \
                         np.roll(vlens,1),np.roll(vlens,2),np.roll(vlens,3)],0)
    localmax = vlens > neighbors
    random_fraction = float(nseeds) / len(angles)
    lengthmask = vlens > 1.5*random_fraction
    vpmask = np.logical_and(localmax,lengthmask)
    vps, vlens, nvps = vps[vpmask,:], vlens[vpmask], np.sum(vpmask)

    # Guided matching to refine the vanishing points
    vcents = np.zeros((nvps,3))
    maxiter = 10
    for i in xrange(nvps):
        vp, vlen = vps[i,:], vlens[i]
        gmtol, gmit, llen = tol, 0, 0
        while vlen!=llen and gmit<maxiter:
            gmit += 1
            linemask = LfromVP(vp,midpts,lineqs,lengths,gmtol)
            vp = VPfromLines(lineqs[linemask,:])
            vlen, llen = np.sum(lengths[linemask]), vlen
        vcent = np.mean(midpts[linemask,:],0)
        vps[i,:], vlens[i], vcents[i,:] = vp, vlen, vcent
    vlens = vlens / total_len

    # eliminate vanishing points without a significant contribution from lines
    keepmask = vlens > 5/90  # "random" line length expected
    vps, vlens, vcents = vps[keepmask,:], vlens[keepmask], vcents[keepmask,:]

    # adjust the sign of vanishing points so that normal associated with vp cross down faces toward camera
    for i in range(len(vlens)):
        vpnorm = np.cross( vps[i,:] , np.dot(tp(Rpr),[0,1,0]) )
        if vpnorm[2] > 0:
            vps[i,:] *= -1
    
    return vps, vlens, vcents, seedlens
        

def VPfromRANSAC(midpts, lineqs, lengths, Rot, maxangle):
    
    # Run a RANSAC loop to determine vanishing points from image lines
    niter, minlen, sinangle = 300, 0.2*np.sum(lengths), np.sin(maxangle*np.pi/180) # RANSAC parameters
    bvps, blens, bcent = np.zeros((0,3)), np.zeros(0), np.zeros((0,3)) # robust estimates
    for i in xrange(niter):
        vp = VPfrom2Lines(lineqs)
        valid, numi, lnumi = True, 0, np.nan
        gmit, mxit = 0, 10
        while valid and numi!=lnumi and gmit<mxit:
            mask = LfromVP(vp,midpts,lineqs,lengths,maxangle)
            vp = VPfromLines(lineqs[mask,:])
            numi, lnumi = np.sum(mask), numi
            valid = validVP(vp, bvps, sinangle, midpts[mask], lengths[mask], minlen, Rot)
            gmit += 1
        vplength = np.sum(lengths[mask]) / np.sum(lengths) # fraction of total length in this vp
        vpcent = np.mean(midpts[mask,:],0)
        if valid: bvps, blens, bcent = np.concatenate( (bvps,[vp]) , 0 ), np.append(blens,vplength), np.concatenate( (bcent,[vpcent]) , 0 )
    sort = np.argsort(blens)[::-1]
    if len(sort) == 0: return bvps, blens, bcent
    bvps, blens, bcent = bvps[sort,:], blens[sort], bcent[sort,:]

    # Remove redundant vanishing points
    idx = 0
    while idx < len(blens):
        vp = bvps[idx,:]
        mask = np.sqrt(1-np.inner(bvps,vp)**2) > sinangle*2
        mask[0] = True
        bvps, blens, bcent = bvps[mask,:], blens[mask], bcent[mask]
        idx += 1
    return bvps, blens, bcent


def VPfrom2Lines(lineqs):
    nlines = lineqs.shape[0]
    i0 = rnd.randint(0,nlines)
    i1 = rnd.randint(0,nlines-1)
    i1 = i1+1 if i1>=i0 else i1
    return geom.normalrows(np.cross(lineqs[i0,:],lineqs[i1,:]))


def VPfromLines(lineqs):
    return geom.smallestSingVector(lineqs)


def validVP(vp, bvps, sinangle, midpts, lengths, minlen, Rot):
    # check to make sure it is not within sinangle of any current vps
    if ( np.sqrt(1-np.inner(bvps,vp)**2) < sinangle/5 ).any(): return False
    # check to make sure none of the contributing lines cross the vp
    #if (geom.vecnorm(geom.vecsub(midpts,vp/vp[2],1))<lengths/2).any(): return False
    # check to make sure the total length contributing isn't too small
    if np.sum(lengths) < minlen: return False
    # check to make sure the incline isn't too great; no more than ~17 degrees
    return np.abs(np.inner(np.dot(Rot,vp),[0,1,0])) < 0.3


def LfromVP(vp,midpts,lineqs,lengths,maxangle):

    # get mask of those within angle tolerance
    line_bear = np.arctan(lineqs[:,0]/lineqs[:,1])
    vp_bear = np.arctan( (midpts[:,1]-vp[1]/vp[2])/(vp[0]/vp[2]-midpts[:,0]) )
    dbear = np.mod(line_bear-vp_bear,np.pi)
    dbear = dbear + (dbear>np.pi/2) * (np.pi-2*dbear)
    anglemask = dbear < maxangle*np.pi/180

    # get mask of those that don't "cross" vanishing point
    vpdist = geom.vecnorm(geom.vecsub(midpts,vp/vp[2],1))
    crossmask = vpdist > lengths/2

    return np.logical_and(anglemask,crossmask)


def LfromLSD(path, img, Kcal):

    # load lines; if not already generated, run LSD
    if not os.path.isdir(os.path.dirname(path)): os.path.mkdir(os.path.dirname(path))
    if not os.path.isfile(path):
        callLSD(path, img)
    lines = loadLines(path)

    # map the line segment endpoints to the image frame
    nlines = lines.shape[0]
    Kinv = alg.inv(Kcal)
    end1 = tp( np.dot( Kinv , np.concatenate( ([lines[:,0]],[lines[:,1]],[np.ones(nlines)]) , 0 ) ) )
    end2 = tp( np.dot( Kinv , np.concatenate( ([lines[:,2]],[lines[:,3]],[np.ones(nlines)]) , 0 ) ) )

    # convert to midpoints, equations, and lengths
    lineqs = np.zeros((nlines,3))
    lineqs[:,0] , lineqs[:,1] = end2[:,1]-end1[:,1] , end1[:,0]-end2[:,0]
    lineqs[:,2] = -np.sum(lineqs*end1,1)
    lineqs = geom.normalrows(lineqs)
    lengths = geom.vecnorm(end1-end2)
    midpts = (end1+end2)/2.0

    # remove lines that are too vertical
    mask = np.abs(lineqs[:,1]/lineqs[:,0]) > np.tan(10*np.pi/180)
    midpts, lineqs, lengths = midpts[mask,:], lineqs[mask,:], lengths[mask]

    return midpts, lineqs, lengths


def loadLines(path):
    data = open(path,'r')
    lines = np.array( [ np.float_(line.strip().split()) for line in data ] )
    data.close()
    return lines


def callLSD(path, img):
    #matlab_path = 'cd(\'/media/DATAPART1/oakland/app/dev-ah/matlab/lsd\'); '
    matlab_path = 'cd(\'/home/jason/Desktop/query/matlab/lsd\'); '
    matlab_lsd  = 'call_lsd(\'' + img + '\',\'' + path + '\'); '
    matlab_call = 'matlab -r \"' + matlab_path + matlab_lsd + 'quit;\"'
    os.system(matlab_call)
    return
