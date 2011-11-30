import geom
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="aaronh"
__date__ ="Nov 29, 2011"

import Image
import time
import numpy as np
from numpy import transpose as tp
import numpy.linalg as alg
import numpy.random as rnd
import render_tags


#####                                                                       #####
#####  Main function for getting qyaw, nyaw from vanishing points analysis  #####
#####                                                                       #####
def getQNyaws(C, qimg, dimg, qsource, nyaw):

    # all VPs within the following degree threshold are considered the same VP
    vp_threshold = 5 # should be similar to the yaw error desired

    # get building face horizontal vp from database image(s) ; compute nyaw
    vp, nyaw = VPNfromDatabase(C, dimg, nyaw, vp_threshold)
    # match vp from query to vp above to compute qyaw
    qyaw = VPQfromQuery(C, qimg, qsource, vp, vp_threshold)

    return qyaw, nyaw


def VPNfromDatabase(C, dimg, nyaw, vp_threshold):

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
    Kl, wRl = viewparam(lsource)
    Kc, wRc = viewparam(csource)
    Kr, wRr = viewparam(rsource)

    # get lines for each database image; image frame equations and segment lengths
    llin, llen = LfromLSD(lpath, limg, Kl)
    clin, clen = LfromLSD(cpath, cimg, Kc)
    rlin, rlen = LfromLSD(rpath, rimg, Kr)

    # get candidate vanishing points from lines
    lvps, lcon = VPfromRANSAC(llin, llen, wRl, vp_threshold)
    cvps, ccon = VPfromRANSAC(clin, clen, wRc, vp_threshold)
    rvps, rcon = VPfromRANSAC(rlin, rlen, wRr, vp_threshold)

    #####  combine candidate vanishing points and nyaw into an estimate of  #####
    #####  the building face's horizontal vanishing point and nyaw          #####

    # increase the confidence of vps from the matched view
    if    view==2 or view==8  : lcon *= 3   # left view
    elif  view==3 or view==9  : ccon *= 2   # less confident of center vps than left/right
    elif  view==4 or view==10 : rcon *= 3   # right view

    # map the vanishing points to the world frame (EDN - east/down/north) and combine all vps
    lvps, cvps, rvps = tp(np.dot(wRl,tp(lvps))), tp(np.dot(wRc,tp(cvps))), tp(np.dot(wRr,tp(rvps)))
    vps, conf = np.concatenate( (lvps,cvps,rvps) , 0 ), np.concatenate((lcon,ccon,rcon))
    nvps = len(conf)

    # get normals and remove vanishing points indicating more than a ~18 degree incline
    normals = np.cross(vps,[0,1,0])
    mask = geom.vecnorm(normals) > 0.95
    vps, normals, conf = vps[mask], geom.normalrows(normals[mask]), conf[mask]
    nvps = len(conf)

    # align vps and remove vps that lead to normals not facing every view
    
    mask = np.bool_(np.zeros(nvps))
    for i in range(nvps):
        if np.dot(tp(wRc),normals[i,:])[2] > 0:
            normals[i,:] *= -1
            vps[i,:] *= -1
        mask[i] = np.dot(tp(wRl),normals[i,:])[2] < 0 and np.dot(tp(wRr),normals[i,:])[2] < 0
    vps, normals, conf = np.compress(mask,vps,0), np.compress(mask,normals,0), np.compress(mask,conf)
    nvps = len(conf)

    # weigh vanishing points by how much they agree with nyaw (if nyaw exists)
    if not np.isnan(nyaw):
        n = np.array([np.sin(nyaw*np.pi/180),0,np.cos(nyaw*np.pi/180)])
        ndots = np.inner(normals,n) # should be close to 1 if they agree
        vpweights = np.float_( ndots > 0.9 ) # less than ~26 degrees away
        conf *= vpweights

    # combine all vanishing points
    bconf, bmask = 0, np.bool_(np.zeros(nvps))
    for i in range(nvps):
        vp = vps[i]
        mask = np.inner(vps,vp) > np.cos(vp_threshold*np.pi/180)
        if np.sum(conf[mask]) > bconf: bconf, bmask = np.sum(conf[mask]), mask
    # weighted sum of inlying vps using the confidence of each
    bvp = geom.normalrows( np.sum( tp(np.tile(conf[bmask],[3,1])) * vps[bmask] , 0 ) )
    bnormal = np.cross(bvp,[0,1,0])
    nyaw = np.arctan2(bnormal[0],bnormal[2])

    return bvp, nyaw


def VPQfromQuery(C, qimg, qsource, vp, vp_threshold):

    # get query vanishing points
    qname = os.path.basename(qimg)
    qpath = os.path.join(C.querydir, 'hires', 'lsd', qname[:-4] + '.lsd')
    Kq, wRq = viewparam(qsource)
    qlin, qlen = LfromLSD(lpath, limg, Kl, wRl)
    qvps, conf = VPfromRANSAC(qlin, qlen, wRq, vp_threshold)

    #####  combine candidate vanishing points and vp from db   #####
    #####  into an estimate of the true query yaw orientation  #####

    # map vanishing points to world frame
    qvps = tp(np.dot(wRq,tp(lvps)))

    # remove vanishing points that deviate in incline from true vp by more than vp_threshold
    qinclines = np.arcsin(np.inner(qvps,[0,1,0])) * 180/np.pi
    incline = np.arcsin(np.inner(vp,[0,1,0])) * 180/np.pi
    mask = np.abs(incline-qinclines) < vp_threshold
    qvps, conf = qvps[mask,:], conf[mask]

    # align vanishing points based on normal
    normals = geom.normalrows(np.cross(qvps,[0,1,0]))
    for i in range(len(conf)):
        if np.dot(tp(wRq),normals[i,:])[2] > 0:
            normals[i,:] *= -1
            qvps[i,:] *= -1
            
    # weigh vanishing points by how close they are to true vp; eliminate if > 30 degrees
    vpsin = np.sqrt( 1 - np.inner(qvps,vp)**2 )
    vpweights = 1 - 2*vpsin if vpsin<0.5 else 0
    conf *= vpweights

    # choose best query vanishing point based on the confidence parameter
    idx = np.argmax(conf)
    qvp, qnorm = qvps[idx,:], normals[idx,:]

    # compute the resulting qyaw by aligning the resulting normal vectors
    cyaw = geom.YPRfromR(wRq)[0] # cell phone yaw
    normal = geom.normalrows(np.cross(vp,[0,1,0]))
    pyaw, qpyaw = np.arctan2(normal[0],normal[2]), np.arctan2(qnorm[0],qnorm[2])
    yawdiff = pyaw - qpyaw
    yawdiff = yawdiff if yawdiff<180 else yawdiff-360
    qyaw = cyaw+yawdiff

    return qyaw


def viewparam(source):
    # camera calibration matrix
    K = geom.cameramat(source.image.size[0], source.image.size[1], source.fov)
    # camera orientation (camera to world)
    R = geom.RfromYPR(source.yaw, source.pitch, source.roll)
    return Kcal, Rot


def VPfromRANSAC(lines, lengths, Rot, vp_threshold):
    
    # Run a RANSAC loop to determine vanishing points from image lines
    niter, minlen, inlerr = 1000, 0.1, np.sqrt(1-np.cos(vp_threshold*np.pi/180)**2) # RANSAC parameters
    bvps, blens = np.zeros((0,3)), np.zeros(0) # robust estimates
    for i in xrange(niter):
        vp = VPfrom2Lines(lines)
        valid, vlen, llen = True, 0, np.nan
        while valid and vlen!=llen:
            mask = np.inner(lines,vp)<inlerr
            vp = VPfromLines(lines[mask])
            vlen, llen = np.sum(lengths[mask]), vlen
            valid = validVP(vp, vlen, bvps, inlerr, minlen, Rot)
        if valid: bvps, blens = np.concatenate( (bvps,[vp]) , 0 ), np.append(blen,vlen)
    return bvps, blens


def VPfrom2Lines(lines):
    nlines = lines.shape[0]
    i0 = rnd.randint(0,nlines)
    i1 = rnd.randint(0,nlines-1)
    i1 = i1+1 if i1>=i0 else i1
    return geom.normalrows(np.cross(lines[i0,:],lines[i1,:])), l


def VPfromLines(lines):
    eigLines = alg.eig( np.dot(tp(lines),lines) )
    return eigLines[1][ : , np.argmax(eigLines[0]) ]


def validVP(vp, vlen, bvps, inlerr, minlen, Rot):
    # check to make sure it is not within inlerr of any current vps
    if ( np.sqrt(1-np.inner(bvps,vp)**2) < inlerr/2 ).any(): return False
    # check to make sure the total length contributing isn't too small
    if vlen < minlen: return False
    # check to make sure the incline isn't too great; no more than ~17 degrees
    return np.abs(np.inner(np.dot(Rot,vp),[0,1,0])) < 0.3


def LfromVP(vp,lines,inlerr):
    return np.inner(lines,vp)<inlerr


def LfromLSD(path, img, Kcal):

    # load lines; if not already generated, run LSD
    if not os.path.isfile(path):
        callLSD(path, img)
    lines = loadLines(path)

    # map the line segment endpoints to the image frame
    nlines = lines.shape[0]
    Kinv = alg.inv(Kcal)
    end1 = tp( np.dot( Kinv , np.concatenate( ([lines[:,1]],[lines[:,0]],[np.ones(nlines)]) , 0 ) ) )
    end2 = tp( np.dot( Kinv , np.concatenate( ([lines[:,3]],[lines[:,2]],[np.ones(nlines)]) , 0 ) ) )

    # convert to equation format and lengths
    lines = np.zeros((nlines,3))
    lines[:,0] , lines[:,1] = end2[:,1]-end1[:,1] , end1[:,0]-end2[:,0]
    lines[:,2] = -np.sum(lines*end1,1)
    lines = geom.normalrows(lines)
    lengths = geom.vecnorm(end1-end2)

    # remove lines that are too vertical
    mask = np.abs(lines[:,1]/lines[:,0]) > 0.2
    lines, lengths = lines[mask], lengths[mask]

    return lines, lengths


def loadLines(path):
    data = open(path,'r')
    lines = np.array( [ np.float_(line.strip().split()) for line in data ] )
    data.close()
    return lines


def callLSD(path, img):
    matlab_path = 'cd(\'/media/DATAPART1/oakland/app/dev-ah/matlab/lsd\'); '
    matlab_lsd  = 'call_lsd(\'' + img + '\',\'' + path + '\'); '
    matlab_call = 'matlab -r \"' + matlab_path + matlab_lsd + 'quit;\"'
    os.system(matlab_call)
    return