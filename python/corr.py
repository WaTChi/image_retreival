#!/usr/bin/env pyhehon

from config import *
import time
import Image, ImageDraw
import render_tags
import numpy as np
import cv
import os

MAX_PIXEL_DEVIATION = 5

def combine_matches(outputFilePaths):
  """Returns dictionary of siftfile => matches"""
  comb = {}
  for res in outputFilePaths:
    detailed = res + '-detailed.npy'
    for image, matches in np.load(detailed):
      assert type(matches) == list
      if image not in comb:
        comb[image] = []
      comb[image].extend(matches)
  return comb

#matches - list of feature match pairs (dict) where each dict {'query':[x,y,scale, rot], 'db':[x,y,scale,rot]}
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
  cv.FindFundamentalMat(pts_q, pts_db, F, status=inliers, param1=MAX_PIXEL_DEVIATION, param2=.99999)
  return F, np.asarray(inliers)[0]

def find_hom(matches):
  pts_q = cv.CreateMat(len(matches), 1, cv.CV_64FC2)
  pts_db = cv.CreateMat(len(matches), 1, cv.CV_64FC2)
  for i, m in enumerate(matches):
    cv.Set2D(pts_q, i, 0, cv.Scalar(*m['query'][:2]))
    cv.Set2D(pts_db, i, 0, cv.Scalar(*m['db'][:2]))
  F = cv.CreateMat(3, 3, cv.CV_64F)
  inliers = cv.CreateMat(1, len(matches), cv.CV_8U)
  cv.SetZero(F)
  cv.SetZero(inliers)
  cv.FindHomography(pts_db, pts_q, F, method=cv.CV_RANSAC, ransacReprojThreshold=MAX_PIXEL_DEVIATION, status=inliers)
  return F, np.asarray(inliers)[0]

def draw_matches(matches, q_img, db_img, out_img, inliers):
  # create image
  assert os.path.exists(q_img)
  assert os.path.exists(db_img)
  a = Image.open(q_img)
  b = Image.open(db_img)
  height = max(a.size[1], b.size[1])
  target = Image.new('RGB', (a.size[0] + b.size[0], height))
  green = np.compress(inliers, matches)
  red = np.compress(map(lambda x: not x, inliers), matches)

  def drawline(match, color):
    db = [match['db'][1] + a.size[0], match['db'][0]]
    draw.line([match['query'][1], match['query'][0]] + db, fill=color, width=3)

  def drawcircle(match):
    draw.ellipse((match['query'][1] - match['query'][2], match['query'][0] - match['query'][2], match['query'][1] + match['query'][2], match['query'][0] + match['query'][2]), outline="hsl(20,100%,50%)")
    draw.ellipse((match['db'][1] + a.size[0] - match['db'][2], match['db'][0] - match['db'][2], match['db'][1] + a.size[0] + match['db'][2], match['db'][0] + match['db'][2]), outline="hsl(20,100%,50%)")

  draw = ImageDraw.Draw(target)

  db = render_tags.TagCollection(os.path.expanduser('~/shiraz/tags.csv'))
  source = render_tags.EarthmineImageInfo(db_img, db_img[:-4] + '.info')
  img = render_tags.TaggedImage(db_img, source, db)
  points = img.map_tags_camera()
  proj_points = []
#  matches = [
#  {'db': (0,10), 'query': (0,20)},
#  {'db': (0,0), 'query': (0,10)},
#  {'db': (10,10), 'query': (10,20)},
#  {'db': (10,0), 'query': (10,10)}
#  ]
  matches = np.compress(inliers, matches)
  H, inliers = find_hom(matches)
  H = np.matrix(np.asarray(H))
  print H
  for (tag, (nulldist, pixel)) in points:
    print pixel
    pixel = H*np.matrix([pixel[0],pixel[1],1]).transpose()
    print pixel
    proj_points.append((tag, (0, pixel)))
  b = img.taggedcopy(points, b)
  print proj_points
  a = img.taggedcopy(proj_points, a)

  target.paste(a, (0,0))
  target.paste(b, (a.size[0],0))
  for x in matches:
    db = x['db']
    querypt = H*np.matrix([db[0], db[1], 1]).transpose()
    q=[abs(querypt[0].item()/querypt[2].item()), abs(querypt[1].item()/querypt[2].item()),0,0]
    print x
    print {'db': x['db'], 'query': q}
    print
    drawline({'db': x['db'], 'query': q}, 'blue')

#  for match in red:
#    drawline(match, 'red')
#    drawcircle(match)
  for match in green:
    drawline(match, 'green')
    drawcircle(match)
  target.save(out_img, 'png')

if __name__ == '__main__':
#  mdir = '/home/ericl/shiraz/'
  mdir = '/home/eric/.gvfs/sftp for ericl on gorgan/home/ericl/shiraz/'
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
  draw_matches(votes, q_img, db_img, '/home/eric/Desktop/out.png', inliers)

# vim: et sw=2
