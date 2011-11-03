# experiments with opencv's support for solving pnp
# doesn't seem to be working very well

import cv
import render_tags
import numpy as np
import numpy.linalg as alg
from numpy import array as arr
from numpy import transpose as tp
from numpy import dot
import geom
import os


# Solves the Perspective N-Point Problem with OpenCV's
# FindExtrinsicCameraParams2().
# returns rotation vector, translation vector as numpy matrices
def solve(C, Q, matches, dbsiftpath, dbimgpath):

  # open EarthMine info
  info = os.path.join(C.infodir, os.path.basename(dbimgpath)[:-4] + '.info')
  em = render_tags.EarthmineImageInfo(dbimgpath, info)
  map3d = C.pixelmap.open(dbsiftpath)

  # find non-None features
  vector = []
  for i, m in enumerate(matches):
    d = m['db']

    # get 3d pt of feature
    feature = map3d[int(d[0]), int(d[1])]
    if not feature:
      continue

    # convert from latlon to meters relative to earthmine camera
    pz, px = geom.lltom(em.lat, em.lon, feature['lat'], feature['lon'])
    py = feature['alt'] - em.alt
    vector.append((m['query'][:2], (px, py, -pz)))

  print vector[0]

  # reference camera matrix
  # f 0 0
  # 0 f 0
  # 0 0 1
  A = cv.CreateMat(3, 3, cv.CV_64F)
  cv.SetZero(A)
  f = 662 # focal len?
  cv.Set2D(A, 0, 0, cv.Scalar(f))
  cv.Set2D(A, 1, 1, cv.Scalar(f))
  cv.Set2D(A, 2, 2, cv.Scalar(1))

  # convert vector to to cvMats
  objectPoints3d = cv.CreateMat(len(vector), 1, cv.CV_64FC3)
  imagePoints2d = cv.CreateMat(len(vector), 1, cv.CV_64FC2)
  for i, (p2d, p3d) in enumerate(vector):
    cv.Set2D(imagePoints2d, i, 0, cv.Scalar(*p2d))
    cv.Set2D(objectPoints3d, i, 0, cv.Scalar(*p3d))

  coeff = cv.CreateMat(4, 1, cv.CV_64F)
  rvec = cv.CreateMat(3, 1, cv.CV_64F)
  tvec = cv.CreateMat(3, 1, cv.CV_64F)
  cv.SetZero(coeff)
  cv.SetZero(rvec)
  cv.SetZero(tvec)

  # since rvec, tvec are zero the initial guess is the earthmine camera
  ret = cv.FindExtrinsicCameraParams2(objectPoints3d, imagePoints2d, A,
    coeff, rvec, tvec, useExtrinsicGuess=False)
  np_rvec = np.matrix(rvec)
  np_tvec = np.matrix(tvec)
  print np_rvec
  print np_tvec
  return np_rvec, np_tvec

def qsolve(C, Q, matches, dbsiftpath, dbimgpath):

    # Get image information.
    info = os.path.join(C.infodir, os.path.basename(dbimgpath)[:-4] + '.info')
    qsource = render_tags.QueryImageInfo(Q.datasource)
    dbsource = render_tags.EarthmineImageInfo(dbimgpath, info)
    map3d = C.pixelmap.open(dbsiftpath)

    # Get 2d pixels and 3d locations of feature inliers
    matches = [(m['query'][0],m['query'][1],m['db'][0],m['db'][1]) for m in matches]
    matches = list(set(matches))
    q2d = [arr([m[0],m[1]]) for m in matches]
    db2d = [arr([m[2],m[3]]) for m in matches]
    db3d = [map3d[int(d[0]),int(d[1])] for d in db2d]
    i = 0
    while i < len(db3d):
        if db3d[i] is None:
            q2d.pop(i)
            db2d.pop(i)
            db3d.pop(i)
        else:
            i = i+1
    olat,olon,oalt = dbsource.lat,dbsource.lon,dbsource.alt
    qlat,qlon = qsource.lat,qsource.lon
    qzx = geom.lltom(olat,olon,qlat,qlon)
    zx = [geom.lltom(olat,olon,d['lat'],d['lon']) for d in db3d]
    y = [dbsource.alt-d['alt'] for d in db3d]
    xyz = [[zx[i][1],y[i],zx[i][0]] for i in range(len(y))]
    print len(xyz)

    # Set K, Rhat
    wx,wy = qsource.pgmsize[0], qsource.pgmsize[1]
    tx,ty = qsource.view_angle[0], qsource.view_angle[1]
    f1, f2 = (wx-1)/(2*np.tan(tx/2)), (wy-1)/(2*np.tan(ty/2))
    f = (f1+f2) / 2
    K = arr([[f,0,(wx-1)/2.0],
             [0,f,(wy-1)/2.0],
             [0,0,1]])
    y,p,r = qsource.yaw, qsource.pitch, qsource.roll
    print [180*y/np.pi,180*p/np.pi,180*r/np.pi]
    Ry = arr([[np.cos(y),0,np.sin(y)],
              [0,1,0],
              [-np.sin(y),0,np.cos(y)]])
    Rx = arr([[1,0,0],
              [0,np.cos(p),-np.sin(p)],
              [0,np.sin(p),np.cos(p)]])
    Rz = arr([[np.cos(r),-np.sin(r),0],
              [np.sin(r),np.cos(r),0],
              [0,0,1]])
    Rhat = dot(Ry,dot(Rx,Rz)) # camera orientation (camera to world)
    tRhat = tp(Rhat)

    # reference camera matrix
    # f 0 cx
    # 0 f cy
    # 0 0 1
    A = cv.CreateMat(3, 3, cv.CV_64F)
    cv.SetZero(A)
    cv.Set2D(A, 0, 0, cv.Scalar(f))
    cv.Set2D(A, 1, 1, cv.Scalar(f))
    cv.Set2D(A, 2, 2, cv.Scalar(1))
    cv.Set2D(A, 0, 2, cv.Scalar((wx-1)/2.0))
    cv.Set2D(A, 1, 2, cv.Scalar((wy-1)/2.0))

    # convert 2d, 3d points to cvMats
    objectPoints3d = cv.CreateMat(len(xyz), 1, cv.CV_64FC3)
    imagePoints2d = cv.CreateMat(len(xyz), 1, cv.CV_64FC2)
    for i in range(len(xyz)):
        cv.Set2D(imagePoints2d, i, 0, cv.Scalar(*q2d[i]))
        cv.Set2D(objectPoints3d, i, 0, cv.Scalar(*xyz[i]))

    # set initial rotation and translation vectors, distortion coefficients
    coeff = cv.CreateMat(4, 1, cv.CV_64F)
    cv.SetZero(coeff)
    rmat = cv.CreateMat(3, 3, cv.CV_64F)
    cv.Set2D(rmat, 0, 0, cv.Scalar(tRhat[0,0]))
    cv.Set2D(rmat, 0, 1, cv.Scalar(tRhat[0,1]))
    cv.Set2D(rmat, 0, 2, cv.Scalar(tRhat[0,2]))
    cv.Set2D(rmat, 1, 0, cv.Scalar(tRhat[1,0]))
    cv.Set2D(rmat, 1, 1, cv.Scalar(tRhat[1,1]))
    cv.Set2D(rmat, 1, 2, cv.Scalar(tRhat[1,2]))
    cv.Set2D(rmat, 2, 0, cv.Scalar(tRhat[2,0]))
    cv.Set2D(rmat, 2, 1, cv.Scalar(tRhat[2,1]))
    cv.Set2D(rmat, 2, 2, cv.Scalar(tRhat[2,2]))
    print 'YPR init for PnP'
    print 180*np.array(geom.ypr_fromR(Rhat))/np.pi
    rvec = cv.CreateMat(3, 1, cv.CV_64F)
    cv.SetZero(rvec)
    cv.Rodrigues2(rmat,rvec) # convert from rotation matrix to Rodrigues vector
    tvec = cv.CreateMat(3, 1, cv.CV_64F)
    cv.SetZero(tvec)

    # solvepnp
    ret = cv.FindExtrinsicCameraParams2(objectPoints3d, imagePoints2d, A,
        coeff, rvec, tvec, useExtrinsicGuess=False)
    np_rvec = np.matrix(rvec).A
    cv.Rodrigues2(rvec,rmat)
    np_rmat = np.transpose(np.matrix(rmat)).A
    np_tvec = np.dot(np_rmat,np.squeeze(np.matrix(tvec).A))
    print 'Rq from PnP'
    print np_rmat
    print 'YPR from PnP'
    print 180*np.array(geom.ypr_fromR(np_rmat))/np.pi

    return np_rmat, np_tvec

def dbsolve(C, Rhat, matches, dbsiftpath, dbimgpath):

    # Get image information.
    info = os.path.join(C.infodir, os.path.basename(dbimgpath)[:-4] + '.info')
    dbsource = render_tags.EarthmineImageInfo(dbimgpath, info)
    map3d = C.pixelmap.open(dbsiftpath)

    # Get 2d pixels and 3d locations of feature inliers
    matches = [(m['query'][0],m['query'][1],m['db'][0],m['db'][1]) for m in matches]
    matches = list(set(matches))
    db2d = [arr([m[2],m[3]]) for m in matches]
    db3d = [map3d[int(d[0]),int(d[1])] for d in db2d]
    i = 0
    while i < len(db3d):
        if db3d[i] is None:
            db2d.pop(i)
            db3d.pop(i)
        else:
            i = i+1
    olat,olon,oalt = dbsource.lat,dbsource.lon,dbsource.alt
    zx = [geom.lltom(olat,olon,d['lat'],d['lon']) for d in db3d]
    y = [dbsource.alt-d['alt'] for d in db3d]
    xyz = [[zx[i][1],y[i],zx[i][0]] for i in range(len(y))]

    # reference camera matrix Kd
    wx,wy = dbsource.image.size
    fov = dbsource.fov
    f = (wx-1)/(2*np.tan(fov/2))
    # f 0 cx
    # 0 f cy
    # 0 0 1
    A = cv.CreateMat(3, 3, cv.CV_64F)
    cv.SetZero(A)
    cv.Set2D(A, 0, 0, cv.Scalar(f))
    cv.Set2D(A, 1, 1, cv.Scalar(f))
    cv.Set2D(A, 2, 2, cv.Scalar(1))
    cv.Set2D(A, 0, 2, cv.Scalar((wx-1)/2.0))
    cv.Set2D(A, 1, 2, cv.Scalar((wy-1)/2.0))

    # convert 2d, 3d points to cvMats
    objectPoints3d = cv.CreateMat(len(xyz), 1, cv.CV_64FC3)
    imagePoints2d = cv.CreateMat(len(xyz), 1, cv.CV_64FC2)
    for i in range(len(xyz)):
        cv.Set2D(imagePoints2d, i, 0, cv.Scalar(*db2d[i]))
        cv.Set2D(objectPoints3d, i, 0, cv.Scalar(*xyz[i]))

    # set initial rotation and translation vectors, distortion coefficients
    coeff = cv.CreateMat(4, 1, cv.CV_64F)
    cv.SetZero(coeff)
    rmat = cv.CreateMat(3, 3, cv.CV_64F)
    tRhat = tp(Rhat)
    cv.Set2D(rmat, 0, 0, cv.Scalar(tRhat[0,0]))
    cv.Set2D(rmat, 0, 1, cv.Scalar(tRhat[0,1]))
    cv.Set2D(rmat, 0, 2, cv.Scalar(tRhat[0,2]))
    cv.Set2D(rmat, 1, 0, cv.Scalar(tRhat[1,0]))
    cv.Set2D(rmat, 1, 1, cv.Scalar(tRhat[1,1]))
    cv.Set2D(rmat, 1, 2, cv.Scalar(tRhat[1,2]))
    cv.Set2D(rmat, 2, 0, cv.Scalar(tRhat[2,0]))
    cv.Set2D(rmat, 2, 1, cv.Scalar(tRhat[2,1]))
    cv.Set2D(rmat, 2, 2, cv.Scalar(tRhat[2,2]))
    rvec = cv.CreateMat(3, 1, cv.CV_64F)
    cv.SetZero(rvec)
    cv.Rodrigues2(rmat,rvec) # convert from rotation matrix to Rodrigues vector
    tvec = cv.CreateMat(3, 1, cv.CV_64F)
    cv.SetZero(tvec)

    #print Rhat
    # solvepnp
    ret = cv.FindExtrinsicCameraParams2(objectPoints3d, imagePoints2d, A,
        coeff, rvec, tvec, useExtrinsicGuess=False)
    np_rvec = np.matrix(rvec).A
    cv.Rodrigues2(rvec,rmat)
    np_rmat = np.transpose(np.matrix(rmat)).A
    np_tvec = np.dot(np_rmat,np.squeeze(np.matrix(tvec).A))

    return np_rmat, np_tvec