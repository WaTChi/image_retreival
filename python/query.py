#!/usr/bin/python
#
# Finds nearest neighbors of a set of features in a database.
# Use Query(...).run() for a single query
# or
# run_parallel(...) to query multiple cells in parallel
#

from config import *
import reader
import util
import corr
import time
import pyflann
import geom
import threading
import numpy as np
import os

PARAMS_DEFAULT = {
  'algorithm': 'kdtree',
  'trees': 1,
  'checks': 1024,
  'dist_threshold': 70000,
# sift or chog or surf
  'descriptor': 'sift',
# euclidean, manhattan, minkowski, hik, hellinger, cs, kl
  'distance_type': 'euclidean',
# use >1 for weighted
  'num_neighbors': 1,
# highest, ransac, ratio, matchonce, filter
  'vote_method': 'filter',
# custom configuration notation
  'confstring': '',
}

class Tree3D:
  def __init__(self, map3d):
    self.map3d = map3d
    self.flann = pyflann.FLANN()
    self.params = PARAMS_DEFAULT.copy()
    self.params['num_neighbors'] = 10
    self.params['checks'] = 4096
    index = '/var/asdf/asdf/asdf' # TODO find place to save index
    if os.path.exists(index):
      INFO('loading index')
      self.flann.load_index(index, map3d)
    else:
      INFO('building index')
      self.flann.build_index(map3d.view(), **self.params)
      INFO('save index?')
      # TODO save index

  def countHigherPtsNear(self, mlat, mlon, malt, threshold):
    dists = []
    querypt = np.ndarray(1, '3float64')
    querypt[0][0] = mlat
    querypt[0][1] = mlon
    querypt[0][2] = malt
    results, dists = self.flann.nn_index(querypt, **self.params)
    results, dists = results[0], dists[0]
    count = 0
    for i, dist in enumerate(dists):
      meters = geom.distance3d6(
        self.map3d[results[i]][0],
        self.map3d[results[i]][1],
        self.map3d[results[i]][2],
        mlat, mlon, malt)
      if meters < 1.0:
        count += 1
      # ignore lower points for rough threshold
      if self.map3d[results[i]][2] < mlat:
        continue
      if meters < threshold:
        count += 1
    return 1 if count > 0 else 0

# I'm not actually sure if the distance function affects
# index compatibility. Someone check please?
def indextype(params):
  dtype = params['distance_type']
  distname = '' if dtype == 'euclidean' else ('-%s' % dtype)
  des = '' if params['descriptor'] == 'sift' else ('-' + params['descriptor'])
  if params['algorithm'] == 'kdtree':
    return 'kdtree%d%s%s' % (params['trees'], distname, des)
  else:
    return '%s%s%s' % (params['algorithm'], distname, des)

def searchtype(params):
  assert len(params) == len(PARAMS_DEFAULT)
  nn = params['num_neighbors']
  nn = '' if nn == 1 else (',nn=%d' % nn)
  vote_method = '' if params['vote_method'] == 'highest' else ',%s' % params['vote_method']
  conf = ''
  if params['confstring']:
    conf = ',%s' % params['confstring']
  return '%s,threshold=%dk,searchparam=%d%s%s%s' % (indextype(params), params['dist_threshold']/1000, params['checks'], nn, vote_method, conf)

def run_parallel(C, Q, cells, outputFilePaths, num_threads=1):
  semaphore = threading.Semaphore(num_threads)
  threads = []
  INFO_TIMING("running queries with up to %d threads" % num_threads)
  for cell, outputFilePath in zip(cells, outputFilePaths):
    thread = Query(C, Q, cell, outputFilePath, semaphore)
    threads.append(thread)
    thread.start()
  for thread in threads:
     thread.join()

# Produces human readable .res file with summary of votes.
# Also produces a .res-detailed.npy file:
#   list of (imagename,
#      list of (vector geom from db,
#              matching vector geom from query))
#   sorted by decreasing number of matching features
#   for plotting, postprocessing, etc.
# Note that the vote method chosen is responsible for this list.
class Query(threading.Thread):
  def __init__(self, C, Q, cell, outfile, barrier=None):
    assert len(C.params) == len(PARAMS_DEFAULT)
    threading.Thread.__init__(self)
    self.qpath = Q.siftpath
    self.cellpath = os.path.join(C.dbdir, cell)
    self.infodir = C.infodir
    self.celldir = C.dbdir
    self.outfile = outfile
    self.params = C.params
    self.barrier = barrier
    self.dump = self.outfile + '-detailed.npy'
    pyflann.set_distance_type(C.params['distance_type'])
    self.reader = reader.get_reader(C.params['descriptor'])

  def run(self):
    if os.path.exists(self.outfile) and os.path.exists(self.dump):
      return
    if self.barrier:
      self.barrier.acquire()
    self.flann = pyflann.FLANN()
    start = time.time()
    dataset, mapping = self._build_index()
    queryset = self.reader.load_file(self.qpath)
    qtime = time.time()
    results, dists = self.flann.nn_index(queryset['vec'], **self.params)
    INFO_TIMING("query took %f seconds" % (time.time() - qtime))
    vtime = time.time()
    sorted_counts = self.vote(queryset, dataset, mapping, results, dists)
    INFO_TIMING("voting took %f seconds" % (time.time() - vtime))
    save_atomic(lambda d: np.save(d, sorted_counts), self.dump)
    votes = [(len(matches), img) for img, matches in sorted_counts]
    def write_votes(d):
      with open(d, 'w') as f:
        for tally in votes:
          f.write("%f\t%s\n" % tally)
    save_atomic(write_votes, self.outfile)
    INFO_TIMING("took %f total" % (time.time() - start))
    if self.barrier:
      self.barrier.release()
    # release memory - in case the Query is still around
    self.flann = None

  def vote(self, queryset, dataset, mapping, results, dists):
    INFO('voting with method %s' % self.params['vote_method'])
    counts = {
      'highest': self._vote_highest,
      'matchonce': self._vote_matchonce,
      'filter': self._vote_filter,
      'ratio': self._vote_spatial_ratio,
      'top_n': self._vote_top_n,
      'ransac': self._vote_ransac,
    }[self.params['vote_method']](queryset, dataset, mapping, results, dists)
    return counts

  def _vote_spatial_ratio(self, queryset, dataset, mapping, results, dists):
    """Like vote_matchonce, but with a spatially aware ratio test.
       Votes must exceed meet ratio with all matches outside of a
       bubble around the original vote location. This requires nn >> 1.
       Params to tweak: ratio threshold, bubble radius"""
    assert self.params['num_neighbors'] > 1
    map3d = self.reader.load_3dmap_for_cell(self.cellpath, dataset, mapping, self.infodir)
    accept, reject, matchonce, ratio = 0, 0, 0, 0
    counts = {} # map from img to counts
    closed = set()
    for i, dist_array in enumerate(dists):
      if results[i][0] in closed:
        matchonce += 1
        reject += 1
        continue
      elif dist_array[0] > self.params['dist_threshold']:
        reject += 1
        continue
      image = mapping[dataset[results[i][0]]['index']]
      coord3d = map3d[results[i][0]]
      passes_ratio_test = False
      if coord3d is None: # no option but to accept
        passes_ratio_test = True
      else:
        dist2 = util.select_dist_exceeds(map3d, mapping, dataset, results[i], coord3d, 10.0, enumerate(dist_array[1:], 1))
        if dist2 is None or dist_array[0]/dist2 < 0.95:
          passes_ratio_test = True
      if not passes_ratio_test:
        ratio += 1
        reject += 1
        continue
      if passes_ratio_test:
        closed.add(results[i][0])
        if image not in counts:
          counts[image] = []
        counts[image].append({'db': dataset[results[i][0]]['geom'].copy(),
                            'query': queryset[i]['geom'].copy()})
        accept += 1
    INFO('accepted %d/%d votes' % (accept, accept + reject))
    INFO('discarded %d vote collisions' % matchonce)
    INFO('discarded %d votes using ratio test' % ratio)
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    return sorted_counts

  def _vote_top_n(self, queryset, dataset, mapping, results, dists):
    """Like vote_matchonce, but up to nn db images per query feature.
       requires more than 1 nearest neighbor for results.
       Note that each image gets 1 vote max."""
    assert self.params['num_neighbors'] > 1
    accept, reject, matchonce = 0, 0, 0
    counts = {} # map from img to counts
    closed = set()
    for i, dist_array in enumerate(dists):
      marked = set()
      for j, dist in enumerate(dist_array):
        if results[i][j] in closed:
          matchonce += 1
          reject += 1
        elif dist < self.params['dist_threshold']:
          closed.add(results[i][j])
          image = mapping[dataset[results[i][j]]['index']]
          if image not in marked:
            if image not in counts:
              counts[image] = []
            counts[image].append({'db': dataset[results[i][j]]['geom'].copy(),
                                'query': queryset[i]['geom'].copy()})
            marked.add(image)
          accept += 1
        else:
          reject += 1
    INFO('accepted %d/%d votes' % (accept, accept + reject))
    if matchonce:
      INFO('discarded %d vote collisions' % matchonce)
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    return sorted_counts

  def false_search(self, queryset):
    self.flann = None # release memory
    # TODO eliminate duplicated build index code
    # TODO eliminate hardcoded path
    falsecellpath = os.path.expanduser('~/shiraz/Research/cells/g=100,r=d=236.6/37.8732916946,-122.279128355')
    falseflann = pyflann.FLANN()
    iname = '%s-%s.%s.index' % (getcellid(falsecellpath), indextype(self.params), np.dtype(self.reader.dtype)['vec'].subdtype[0].name)
    index = getfile(falsecellpath, iname)
    dataset, mapping = self.reader.load_cell(falsecellpath)
    if os.path.exists(index):
      falseflann.load_index(index, dataset['vec'])
    else:
      falseflann.build_index(dataset['vec'], **self.params)
      for out in getdests(falsecellpath, iname):
        save_atomic(lambda d: falseflann.save_index(d), out)
    dists = []
    r, dists = falseflann.nn_index(queryset['vec'], **self.params)
    return r, dists

  def _vote_filter(self, queryset, dataset, mapping, results, dists):
    """Votes must beat false votes in another cell."""
    assert self.params['num_neighbors'] == 1
    counts = {} # map from img to counts
    closed = set()
    closed2 = set()
    accept, reject, matchonce, vs, c2 = 0, 0, 0, 0, 0
    results2, contest = self.false_search(queryset)
    for i, dist in enumerate(dists):
      if dist > self.params['dist_threshold']:
        reject += 1
      elif dist > contest[i] and results2[i] in closed2:
        reject += 1
        c2 += 1
      elif dist > contest[i]:
        closed2.add(results2[i])
        reject += 1
        vs += 1
      elif results[i] in closed:
        reject += 1
        matchonce += 1
      else:
        closed.add(results[i])
        accept += 1
        img = mapping[dataset[results[i]]['index']]
        if img not in counts:
          counts[img] = []
        counts[img].append({'db': dataset[results[i]]['geom'].copy(),
                            'query': queryset[i]['geom'].copy()})
    INFO('accepted %d/%d votes' % (accept, accept + reject))
    if matchonce:
      INFO('discarded %d vote collisions' % matchonce)
    if vs:
      INFO('discarded %d losing votes' % vs)
    if c2:
      INFO('%d votes escaped filtering' % c2)
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    return sorted_counts

  def _vote_matchonce(self, queryset, dataset, mapping, results, dists):
    """Like vote highest, but each db feature is matchonceed to 1 match"""
    assert self.params['num_neighbors'] == 1
    counts = {} # map from img to counts
    closed = set()
    accept, reject, matchonce = 0, 0, 0
    for i, dist in enumerate(dists):
      if dist > self.params['dist_threshold']:
        reject += 1
      elif results[i] in closed:
        reject += 1
        matchonce += 1
      else:
        closed.add(results[i])
        accept += 1
        img = mapping[dataset[results[i]]['index']]
        if img not in counts:
          counts[img] = []
        counts[img].append({'db': dataset[results[i]]['geom'].copy(),
                            'query': queryset[i]['geom'].copy()})
    INFO('accepted %d/%d votes' % (accept, accept + reject))
    if matchonce:
      INFO('discarded %d vote collisions' % matchonce)
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    return sorted_counts

  def _vote_highest(self, queryset, dataset, mapping, results, dists):
    assert self.params['num_neighbors'] == 1
    counts = {} # map from img to counts
    accept, reject = 0, 0
    for i, dist in enumerate(dists):
      if dist < self.params['dist_threshold']:
        accept += 1
        img = mapping[dataset[results[i]]['index']]
        if img not in counts:
          counts[img] = []
        counts[img].append({'db': dataset[results[i]]['geom'].copy(),
                            'query': queryset[i]['geom'].copy()})
      else:
        reject += 1
    INFO('accepted %d/%d votes' % (accept, accept + reject))
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    return sorted_counts

  def _vote_ransac(self, queryset, dataset, mapping, results, dists):
    if self.params['num_neighbors'] > 1:
      sorted_counts = self._vote_top_n(queryset, dataset, mapping, results, dists)
    else:
      sorted_counts = self._vote_matchonce(queryset, dataset, mapping, results, dists)
    # filters out outliers from counts until
    # filtered_votes(ith image) > votes(jth image) for all j != i
    # and returns top 10 filtered results
    filtered = {}
    bound = -1
    num_filt = 0
    for siftfile, matches in sorted_counts:
      if len(matches) < bound or num_filt > 10:
        INFO('stopped after filtering %d' % num_filt)
        break
      num_filt += 1
      F, inliers = corr.find_corr(matches)
      bound = max(sum(inliers), bound)
      pts = np.ndarray(len(matches), np.object)
      pts[0:len(matches)] = matches
      if any(inliers):
        filtered[siftfile] = list(np.compress(inliers, pts))
    rsorted_counts = sorted(filtered.iteritems(), key=lambda x: len(x[1]), reverse=True)
    if not rsorted_counts:
      INFO('W: ransac rejected everything, not filtering')
      return sorted_counts
    return rsorted_counts

  def _build_index(self):
    start = time.time()
    iname = '%s-%s.%s.index' % (getcellid(self.cellpath), indextype(self.params), np.dtype(self.reader.dtype)['vec'].subdtype[0].name)
    index = getfile(self.cellpath, iname)
    dataset, mapping = self.reader.load_cell(self.cellpath)
    INFO_TIMING("dataset load took %f seconds" % (time.time() - start))
    if os.path.exists(index):
      s = time.time()
      self.flann.load_index(index, dataset['vec'])
      INFO_TIMING("index load took %f seconds" % (time.time() - s))
      return dataset, mapping
    INFO('creating %s' % iname)
    start = time.time()
    INFO(self.flann.build_index(dataset['vec'], **self.params))
    INFO_TIMING("index creation took %f seconds" % (time.time() - start))
    for out in getdests(self.cellpath, iname):
      save_atomic(lambda d: self.flann.save_index(d), out)
    return dataset, mapping

# vim: et sw=2
