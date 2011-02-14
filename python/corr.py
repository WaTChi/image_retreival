#!/usr/bin/env python

from config import *
import random
import time
import Image, ImageDraw
import render_tags
import numpy as np
import cv
import os

MAX_PIXEL_DEVIATION = 5
CONFIDENCE_LEVEL = .99999
ROT_THRESHOLD_RADIANS = 0.2 # .1 ~ 5 deg
best_rot = 0

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

MAX_SPATIAL_ERROR = 0
def getSpatiallyOrdered(matches, axis, inliers):
  """return length of max increasing subseq"""
  indices = map(lambda (i,m): (m['db'][axis], m['query'][axis], i), filter(lambda (j,m): inliers[j], enumerate(matches)))
  if not indices:
    return inliers
  indices.sort() # by db order
  # now find subseq by query order
  L = [None]*len(indices)
  def edgesof(j):
    for i in range(j):
      if indices[i][1] <= indices[j][1] + MAX_SPATIAL_ERROR:
        yield i
  def merge((c1, l1), (c2, l2)):
    return (c1 + c2, l1 + l2)
  for j in range(len(indices)):
    E = list(edgesof(j))
    if E:
      L[j] = merge((1, [indices[j][2]]), max([L[i] for i in E]))
    else:
      L[j] = (1, [indices[j][2]])
  best = set(max(L)[1])
  for i in range(len(inliers)):
    if i not in best:
      inliers[i] = False
  return inliers

def rot_delta(m, correction=0):
  a = m['query'][3] + correction
  b = m['db'][3]
  rot = 2*np.pi
  return min((a-b) % rot, (b-a) % rot)

#matches - list of feature match pairs (dict) where each dict {'query':[x,y,scale, rot], 'db':[x,y,scale,rot]}
def find_corr(matches, hom=False):
  matches = list(matches)
  F = cv.CreateMat(3, 3, cv.CV_64F)
  cv.SetZero(F)
  if not matches or (hom and len(matches) < 4):
    return F, []
  inliers = cv.CreateMat(1, len(matches), cv.CV_8U)
  cv.SetZero(inliers)
  pts_q = cv.CreateMat(len(matches), 1, cv.CV_64FC2)
  pts_db = cv.CreateMat(len(matches), 1, cv.CV_64FC2)
  for i, m in enumerate(matches):
    cv.Set2D(pts_q, i, 0, cv.Scalar(*m['query'][:2]))
    cv.Set2D(pts_db, i, 0, cv.Scalar(*m['db'][:2]))

# ransac for fundamental matrix. rotation filtering
# TODO multiple RANSAC to get smaller/larger features
  if not hom:
    cv.FindFundamentalMat(pts_q, pts_db, F, status=inliers, param1=MAX_PIXEL_DEVIATION, param2=CONFIDENCE_LEVEL)
    inliers = np.asarray(inliers)[0]
    global best_rot
# TODO use homography to find correct orientation
# this will fix assumption that db,query have same roll
    best_rot = (-9999, 0)
    for i, m in enumerate(matches):
      if inliers[i]:
        if abs(rot_delta(m, best_rot[1])) > ROT_THRESHOLD_RADIANS:
          inliers[i] = False
## reduces performance for GT %
## TODO check if it increases localization %
#    inliers = getSpatiallyOrdered(matches, 0, inliers)
#    inliers = getSpatiallyOrdered(matches, 1, inliers)
    return F, inliers

  # homography only. no rotation check
  cv.FindHomography(pts_db, pts_q, F, method=cv.CV_RANSAC, ransacReprojThreshold=MAX_PIXEL_DEVIATION, status=inliers)
  return F, np.asarray(inliers)[0]

def draw_matches(matches, q_img, db_img, out_img, inliers, showLine=True, showtag=True):
  # create image
  assert os.path.exists(q_img)
  assert os.path.exists(db_img)
  a = Image.open(q_img)
  if a.mode != 'RGB':
    a = a.convert('RGB')
  scale = 1
  portrait = a.size[0] < a.size[1]
  if a.size[0] > 768 or a.size[1] > 512:
    newy = 512
    scale = float(newy)/a.size[1]
    newx = a.size[0]*scale
    INFO('resizing image %s => %s' % (str(a.size), str((newx,newy))))
    a = a.resize((newx,newy), Image.ANTIALIAS)
    # XXX TODO rework rescale
    if portrait:
      scale = float(newy)/840
#      scale = 1
  assert a.mode == 'RGB'
  b = Image.open(db_img)
  height = max(a.size[1], b.size[1])
  target = Image.new('RGBA', (a.size[0] + b.size[0], height))

  def drawline(match, color='hsl(20,100%,50%)', w=3):
    db = [match['db'][1] + a.size[0], match['db'][0]]
    draw.line([match['query'][1], match['query'][0]] + db, fill=color, width=w)

  def drawcircle(match, col='hsl(20,100%,50%)'):
    draw.ellipse((match['query'][1] - match['query'][2], match['query'][0] - match['query'][2], match['query'][1] + match['query'][2], match['query'][0] + match['query'][2]), outline=col)
    draw.ellipse((match['db'][1] + a.size[0] - match['db'][2], match['db'][0] - match['db'][2], match['db'][1] + a.size[0] + match['db'][2], match['db'][0] + match['db'][2]), outline=col)

  draw = ImageDraw.Draw(target)

  db = render_tags.TagCollection(os.path.expanduser('~/shiraz/Research/app/dev/tags.csv'))
  source = render_tags.EarthmineImageInfo(db_img, db_img[:-4] + '.info')
  img = render_tags.TaggedImage(db_img, source, db)
  points = img.map_tags_camera()
  proj_points = []
  oldmatches = matches
  oldinliers = inliers
  matches = np.compress(inliers, matches)
  H, inliers = find_corr(matches, hom=True)
  H = np.matrix(np.asarray(H))
  tagmatches = []

  #TODO: include outliers from homography in red, refactor stuff. the oldinlier/inlier stuff is confusing
  green = np.compress(inliers, matches).tolist()
  red = np.compress(map(lambda x: not x, oldinliers), oldmatches).tolist()
  red += np.compress(map(lambda x: not x, inliers), matches).tolist()

  for g in green:
      g['query'][0]*=scale
      g['query'][1]*=scale

  for g in red:
      g['query'][0]*=scale
      g['query'][1]*=scale

  # confusing geometry. x and y switch between the reprs.
  for (tag, (nulldist, pixel)) in points:
    x = pixel[1]
    y = pixel[0]
    dest = H*np.matrix([x,y,1]).transpose()
    try:
      dest = tuple(map(int, (dest[0].item()/dest[2].item(), dest[1].item()/dest[2].item())))
    except ZeroDivisionError:
      dest = (0,0)
    except ValueError:
      dest = (0,0)
    tagmatches.append({'db': [x, y, 10], 'query': [dest[0], dest[1], 10]})
    dest = (dest[1]*scale, dest[0]*scale)
    proj_points.append((tag, (0, dest)))

  target.paste(a, (0,0))
  target.paste(b, (a.size[0],0))

  def colorize(theta):
    if theta < .1:
      return 'green'
    elif theta < .2:
      return 'blue'
    elif theta < .5:
      return 'purple'
    elif theta < 1.5:
      return 'orange'
    else:
      return 'red'

  if showLine:
      for match in red:
        drawline(match, 'red', w=1)
        drawcircle(match, colorize(rot_delta(match, best_rot[1])))
      for match in green:
        drawline(match, 'green', w=2)
        drawcircle(match, colorize(rot_delta(match, best_rot[1])))

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
  if showtag:
    target.paste(tagfilled, mask=tags)

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
