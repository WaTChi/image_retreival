#!/usr/bin/env python

from config import *
import time
import numpy as np
import cv
import os

CONFIDENCE_LEVEL = 0.99
MAX_PIXEL_DEVIATION = 3

def find_corr(matches):
  pts_q = cv.CreateMat(len(matches), 1, cv.CV_64FC2)
  pts_db = cv.CreateMat(len(matches), 1, cv.CV_64FC2)
  for i, m in enumerate(matches):
    cv.Set2D(pts_q, i, 0, cv.Scalar(*m['query'][:2]))
    cv.Set2D(pts_db, i, 0, cv.Scalar(*m['db'][:2]))
  F = cv.CreateMat(3, 3, cv.CV_64F)
  inliers = cv.CreateMat(1, len(matches), cv.CV_8U)
  cv.SetZero(F)
  cv.SetZero(inliers)
  ret = cv.FindFundamentalMat(pts_q, pts_db, F, status=inliers, param1=MAX_PIXEL_DEVIATION, param2=CONFIDENCE_LEVEL)
#  if not ret:
#    raise Exception, "No correspondence found."
  return F, np.asarray(inliers)[0]

def draw_matches(detailed_resfile, q_img, db_img, out_img, inliers):
  assert os.path.exists(q_img)
  assert os.path.exists(db_img)
  pass

if __name__ == '__main__':
  mdir = '/home/ericl/shiraz/'
  res = mdir + 'Research/results/query3/matchescells(g=100,r=d=236.6),query3,kdtree1,threshold=70k,searchparam=1024/DSC_8848,37.872682,-122.268556sift.txt,37.8714488812,-122.266998471,193.626818312.res-detailed.npy'
  votes = np.load(res)[0][1] # top image
  start = time.time()
  F, inliers = find_corr(votes)
  INFO('F =\n'+str(np.asarray(F)))
  n = sum(inliers)
  d = len(inliers)
  INFO("RANSAC left %d/%d (%d%%) inliers" % (n,d,100*n/d))
  INFO_TIMING('RANSAC took %f seconds' % (time.time() - start))
  q_img = mdir + 'query3/DSC_8848,37.872682,-122.268556.pgm'
  db_img = mdir + 'Research/collected_images/earthmine-fa10.1/37.871955,-122.270829/37.8726500925,-122.268473956-0010.jpg'
  draw_matches(res, q_img, db_img, '/home/ericl/out.png', inliers)

# vim: et sw=2
