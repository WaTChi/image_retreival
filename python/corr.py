#!/usr/bin/env pyhehon

from config import *
import random
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
  q_img = q_img.replace('.pgm', '.jpg')
  a = Image.open(q_img)
  if a.mode != 'RGB':
    a = a.convert('RGB')
  if a.size > (768, 512):
    INFO('resizing image %s => %s' % (str(a.size), '(768, 512)'))
    a = a.resize((768, 512), Image.ANTIALIAS)
  assert a.mode == 'RGB'
  b = Image.open(db_img)
  height = max(a.size[1], b.size[1])
  target = Image.new('RGBA', (a.size[0] + b.size[0], height))

  def drawline(match, color='hsl(20,100%,50%)', w=1):
    db = [match['db'][1] + a.size[0], match['db'][0]]
    draw.line([match['query'][1], match['query'][0]] + db, fill=color, width=w)

  def drawcircle(match, col='hsl(20,100%,50%)'):
    draw.ellipse((match['query'][1] - match['query'][2], match['query'][0] - match['query'][2], match['query'][1] + match['query'][2], match['query'][0] + match['query'][2]), outline=col)
    draw.ellipse((match['db'][1] + a.size[0] - match['db'][2], match['db'][0] - match['db'][2], match['db'][1] + a.size[0] + match['db'][2], match['db'][0] + match['db'][2]), outline=col)

  draw = ImageDraw.Draw(target)

  db = render_tags.TagCollection(os.path.expanduser('../tags.csv'))
  source = render_tags.EarthmineImageInfo(db_img, db_img[:-4] + '.info')
  img = render_tags.TaggedImage(db_img, source, db)
  points = img.map_tags_camera()
  proj_points = []
  oldmatches = matches
  oldinliers = inliers
  matches = np.compress(inliers, matches)
  H, inliers = find_hom(matches)
  H = np.matrix(np.asarray(H))
  tagmatches = []

  green = np.compress(inliers, matches)
  red = np.compress(map(lambda x: not x, oldinliers), oldmatches)

  # deeply confusing geometry. x and y switch between the reprs.
  for (tag, (nulldist, pixel)) in points:
    x = pixel[1]
    y = pixel[0]
    dest = H*np.matrix([x,y,1]).transpose()
    dest = tuple(map(int, (dest[0].item()/dest[2].item(), dest[1].item()/dest[2].item())))
    tagmatches.append({'db': [x, y, 10], 'query': [dest[0], dest[1], 10]})
    dest = (dest[1], dest[0])
    proj_points.append((tag, (0, dest)))

  target.paste(a, (0,0))
  target.paste(b, (a.size[0],0))

  for match in red:
    drawline(match)
    drawcircle(match)
  for match in green:
    drawline(match, 'yellow')
    drawcircle(match, 'yellow')
  # ImageDraw :(
  a2 = img.taggedcopy(proj_points, a)
  b2 = img.taggedcopy(points, b)
  a = Image.new('RGBA', (a.size[0], height))
  b = Image.new('RGBA', (b.size[0], height))
  a = img.taggedcopy(proj_points, a)
  b = img.taggedcopy(points, b)
  tags = Image.new('RGBA', (a.size[0] + b.size[0], height))
  tagfilled = Image.new('RGBA', (a.size[0] + b.size[0], height))
  tags.paste(a, (0,0))
  tags.paste(b, (a.size[0],0))
  tagfilled.paste(a2, (0,0))
  tagfilled.paste(b2, (a.size[0],0))
  target.paste(tagfilled, mask=tags)

#  for match in tagmatches:
#    rand = 'hsl(%d,100%%,50%%)' % (100 + int(155*random.random()))
#    drawline(match, rand)
#    drawcircle(match)
#
#  for match in oldmatches:
#    dest = H*np.matrix([match['db'][0],match['db'][1],1]).transpose()
#    dest = tuple(map(int, (dest[0].item()/dest[2].item(), dest[1].item()/dest[2].item())) + [10])
#    match['query'] = dest
#    drawline(match, 'red')
#    drawcircle(match)
  target.save(out_img, 'jpeg', quality=90)

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
