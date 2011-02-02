#!/usr/bin/env python

import cv
import earthMine as em
import pixels

img = "37.8695551919,-122.266734533-0004sift.txt"
v = {'view-location':{'lat':37.86955519192, 'lon':-122.26673453337, 'alt': 62.466999053955078}}

# read in points
px = pixels.PixelMap('/home/zhangz/shiraz/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829')
rawlocs = px.open(img)
#filter out ones w/o 3d points
locsar = filter(lambda x: x[1], rawlocs.items())

# change 3d coordinate systems
c = map(lambda x: ({'x':x[0][0], 'y':x[0][1]}, x[1]), locsar)
locpairs = em.ddImageLocstoLPT(v, c)

pts2d = map(lambda x: x[0], locpairs)
pts3d = map(lambda x: tuple(x[1]), locpairs)
positobj = cv.CreatePOSITObject(pts3d)

rotMat, transVec = cv.POSIT(positobj, pts2d, 0.008, (cv.CV_TERMCRIT_EPS, 0, 0.01))

#change 3d coordinate data format
pts3d_mat = cv.CreateMat(len(pts3d), 1, cv.CV_64FC3)
for i, m in enumerate(pts3d):
    cv.Set2D(pts3d_mat, i, 0, cv.Scalar(*m))

#change rotMat to cvMat
rotMatCV = cv.CreateMat(3, 3, cv.CV_64F)
for i in range(0, 3):
    for j in range(0, 3):
        cv.Set2D(rotMatCV,i, j, rotMat[i][j])

#change transVec to cvMat
transVecCV = cv.CreateMat(1, 3, cv.CV_64F)
for i in range(0, 3):
    transVecCV[0, i] = transVec[i]

#camera matrix
cameratrans = cv.CreateMat(3, 3, cv.CV_64F)
cv.SetIdentity(cameratrans)


#convert rotMatrix to rotVector
rotVec = cv.CreateMat(3, 1, cv.CV_64F)
cv.Rodrigues2(rotMatCV, rotVec)

#distCoefs
distCoef = cv.CreateMat(4, 1, cv.CV_64F)
cv.SetZero(distCoef)

d2 = cv.CreateMat(pts3d_mat.rows, 1, cv.CV_64FC2)
cv.ProjectPoints2(pts3d_mat, rotVec, transVecCV, cameratrans, distCoef, d2)
for i in range(700, d2.rows):
    print d2[i,0]

# vim: et sw=2
