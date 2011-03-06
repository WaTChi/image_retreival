#!/usr/bin/env python

from config import *
import random
import pyflann
import time
import Image, ImageDraw
import render_tags
import numpy as np
import cv
import os

fail, ok = 0, 0
MAX_PIXEL_DEVIATION = 5
FALLBACK_PIXEL_DEVIATIONS = [2,1]
CONFIDENCE_LEVEL = .9999999
BAD_HOMOGRAPHY_DET_THRESHOLD = .005
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

def hashmatch(m):
  d = m['db']
  q = m['query']
  o = 0
  for i in d:
    o += i
  for i in q:
    o += i
  return o

# eqv of 70k, euclidean, matchonce
def rematch(reader, querysift, dbsift):
#  INFO('Attempting rematch between %s and %s' % (querysift, dbsift))
  start = time.time()
  q = reader.load_file(querysift)
  db = reader.load_file(dbsift)
  flann = pyflann.FLANN()
  results, dists = flann.nn(db['vec'], q['vec'], 1, algorithm='linear')
  matches = []
  closed = set()
  for i in range(0,len(results)):
    if dists[i] < 70000 and results[i] not in closed:
      closed.add(results[i])
      atom = {'db': db[results[i]]['geom'].copy(),
                      'query': q[i]['geom'].copy()}
      matches.append(atom)
  INFO_TIMING("rematch took %f" % (time.time() - start))
  return matches

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
# if you use ransac_pass, the inliers returned will come with
#   a reduced set of matches
# success - array len > 0 : success[0] indicates hom. estimated success
# undefined behavior in ransac case
def find_corr(matches, hom=False, ransac_pass=True, success=[False]):
  if hom:
    if ransac_pass:
      F, inliers = _find_corr(matches, MAX_PIXEL_DEVIATION=50, rotation_filter_only=False, ROT_THRESHOLD_RADIANS=30*np.pi/180)
      matches = np.compress(inliers, matches)
    return (matches,) + _find_corr(matches, hom=True, success=success, MAX_PIXEL_DEVIATION=35, FALLBACK_PIXEL_DEVIATIONS=[15, 9, 5, 3, 1])
  else:
    return _find_corr(matches, success=success)

def _find_corr(matches, hom=False, success=[False], MAX_PIXEL_DEVIATION=MAX_PIXEL_DEVIATION, FALLBACK_PIXEL_DEVIATIONS=FALLBACK_PIXEL_DEVIATIONS, rotation_filter_only=False, ROT_THRESHOLD_RADIANS=ROT_THRESHOLD_RADIANS):
  global fail, ok
  matches = list(matches)
  F = cv.CreateMat(3, 3, cv.CV_64F)
  cv.SetZero(F)
  if not matches or (hom and len(matches) < 4):
    success[0] = False
    fail += 1
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
    if rotation_filter_only:
      inliers = [1 for i in inliers]
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

### try rounds of homography calculations ###
  i = 0
  while i < len(FALLBACK_PIXEL_DEVIATIONS):
    if not isHomographyGood(F):
      cv.FindHomography(pts_db, pts_q, F, method=cv.CV_RANSAC, ransacReprojThreshold=FALLBACK_PIXEL_DEVIATIONS[i], status=inliers)
      i += 1
    else:
#      if i > 0:
#        INFO('SUCCESS in revising homography')
      break
  if i >= len(FALLBACK_PIXEL_DEVIATIONS):
#    INFO('FAILED to revise homography')
    cv.FindHomography(pts_db, pts_q, F, method=cv.CV_LMEDS, status=inliers)
    if isHomographyGood(F):
#      INFO('(!!!) LMEDS worked where RANSAC did not')
      success[0] = True
      ok += 1
    else:
      fail += 1
      success[0] = False
#      INFO('LMEDS also failed')
  else:
    success[0] = True
    ok += 1
  INFO('homography success %d/%d' % (ok,(ok+fail)))

  return F, np.asarray(inliers)[0]

def isHomographyGood(H):
  pts = [(384,256), (384,361), (489,256)]
  det = np.linalg.det(H)
  if det < BAD_HOMOGRAPHY_DET_THRESHOLD:
    return False
  scale = 1
  dests = []
  for (y,x) in pts:
    dest = H*np.matrix([x,y,1]).transpose()
    try:
      dest = tuple(map(int, (dest[0].item()/dest[2].item(), dest[1].item()/dest[2].item())))
    except ZeroDivisionError:
      dest = (0,0)
    except ValueError:
      dest = (0,0)
    dest = (dest[1]*scale, dest[0]*scale)
    dests.append(dest)
  v1len = np.sqrt((dests[0][0] - dests[1][0])**2 + (dests[0][1] - dests[1][1])**2)
  v2len = np.sqrt((dests[0][0] - dests[2][0])**2 + (dests[0][1] - dests[2][1])**2)
  length_ratio = v1len / v2len
  if length_ratio > 5 or 1/length_ratio > 5:
    return False
  a = (dests[1][0] - dests[0][0], dests[1][1] - dests[0][1])
  b = (dests[2][0] - dests[0][0], dests[2][1] - dests[0][1])
  aLen = (a[0]**2 + a[1]**2)**.5
  bLen = (b[0]**2 + b[1]**2)**.5
  angle_delta = np.arccos(np.dot(a, b)/(aLen*bLen))*180.0/np.pi
  if angle_delta < 30 or angle_delta > 160:
    return False
  return True

# returns H, inliers
def draw_matches(matches, q_img, db_img, out_img, showLine=True, showtag=True, showHom=False, success=[False]):
  # create image
  INFO(q_img)
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
    INFO_TIMING('resizing image %s => %s' % (str(a.size), str((newx,newy))))
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

  def xdrawline((start,stop), color='hsl(20,100%,50%)', off=0):
    start = [start[0] + off, start[1]]
    stop = [stop[0] + off, stop[1]]
    draw.line(start + stop, fill=color, width=8)

  def xdrawcircle((y,x), col='hsl(20,100%,50%)', off=0):
    r = 50
    draw.ellipse((y-r+off, x-r, y+r+off, x+r), outline=col)

  def drawcircle(match, col='hsl(20,100%,50%)'):
    draw.ellipse((match['query'][1] - match['query'][2], match['query'][0] - match['query'][2], match['query'][1] + match['query'][2], match['query'][0] + match['query'][2]), outline=col)
    draw.ellipse((match['db'][1] + a.size[0] - match['db'][2], match['db'][0] - match['db'][2], match['db'][1] + a.size[0] + match['db'][2], match['db'][0] + match['db'][2]), outline=col)

  draw = ImageDraw.Draw(target)

  db = render_tags.TagCollection(os.path.expanduser('/media/DATAPART2/Research/app/dev/tags.csv'))
  source = render_tags.EarthmineImageInfo(db_img, db_img[:-4] + '.info')
  img = render_tags.TaggedImage(db_img, source, db)
  points = img.map_tags_camera()
  proj_points = []
  rsc_matches, H, inliers = find_corr(matches, hom=True, ransac_pass=True, success=success)
  H = np.matrix(np.asarray(H))
  tagmatches = []

  green = np.compress(inliers, rsc_matches).tolist()
  notred = set([hashmatch(m) for m in green])
  red = []
  for m in matches:
    if hashmatch(m) not in notred:
      red.append(m)

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
#      for match in red:
#        drawline(match, 'red', w=1)
#        drawcircle(match, colorize(rot_delta(match, best_rot[1])))
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

  if showHom:
    pts = [(384,256), (384,361), (489, 256)]
    dests = []
    for (y,x) in pts:
      dest = H*np.matrix([x,y,1]).transpose()
      try:
        dest = tuple(map(int, (dest[0].item()/dest[2].item(), dest[1].item()/dest[2].item())))
      except ZeroDivisionError:
        dest = (0,0)
      except ValueError:
        dest = (0,0)
      dest = (dest[1]*scale, dest[0]*scale)
      dests.append(dest)

    while dests[0][0] < 100:
      dests[0] = dests[0][0] + 100, dests[0][1]
      dests[1] = dests[1][0] + 100, dests[1][1]
      dests[2] = dests[2][0] + 100, dests[2][1]
    while dests[0][0] > a.size[0] - 100:
      dests[0] = dests[0][0] - 100, dests[0][1]
      dests[1] = dests[1][0] - 100, dests[1][1]
      dests[2] = dests[2][0] - 100, dests[2][1]
    while dests[0][1] > a.size[1] - 100:
      dests[0] = dests[0][0], dests[0][1] - 100
      dests[1] = dests[1][0], dests[1][1] - 100
      dests[2] = dests[2][0], dests[2][1] - 100
    while dests[0][1] < 100:
      dests[0] = dests[0][0], dests[0][1] + 100
      dests[1] = dests[1][0], dests[1][1] + 100
      dests[2] = dests[2][0], dests[2][1] + 100

    xdrawline((pts[0], pts[1]), 'violet', off=a.size[0])
    xdrawline((pts[0], pts[2]), 'yellow', off=a.size[0])
    xdrawline((dests[0], dests[1]), 'violet')
    xdrawline((dests[0], dests[2]), 'yellow')
    xdrawcircle(dests[0], 'yellow')
    xdrawcircle(pts[0], 'yellow', off=a.size[0])
    det = np.linalg.det(H)
    detstr = " det H = %s" % det
    truematch = isHomographyGood(H)
    guess = " guess = %s match" % truematch
    draw.rectangle([(0,0), (150,20)], fill='black')
    draw.text((0,0), detstr)
    draw.text((0,10), guess)

  target.save(out_img, 'jpeg', quality=90)
  return H, inliers

# vim: et sw=2
