#!/usr/bin/env python

import pixels
import earthMine as em
import cv

img = "37.8695551919,-122.266734533-0004sift.txt"
v = {'view-location':{'lat':37.86955519192, 'lon':-122.26673453337, 'alt': 62.466999053955078}}

# read in points
px = pixels.PixelMap('/home/zhangz/shiraz/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829')
rawlocs = px.open(img)
locsar = filter(lambda x: x[1], rawlocs.items())

# change 3d coordinate systems
c = map(lambda x: ({'x':x[0][0], 'y':x[0][1]}, x[1]), locsar)
locpairs = em.ddImageLocstoLPT(v,c)

pts2d = map(lambda x: x[0][:2], locpairs)
pts3d = map(lambda x: tuple(x[1]), locpairs)
positobj = cv.CreatePOSITObject(pts3d)

rotationMatrix, translationVector = \
  cv.POSIT(positobj, pts2d, 20, (cv.CV_TERMCRIT_EPS, 0, 0.01))

print rotationMatrix
print translationVector

# vim: et sw=2
