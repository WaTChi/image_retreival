import cv
import render_tags
import numpy as np
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

# vim: et sw=2
