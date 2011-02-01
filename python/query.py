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
from multiprocessing import cpu_count
import time
import pyflann
import threading
import numpy as np
import os

PARAMS_DEFAULT = {
  'algorithm': 'kdtree',
  'trees': 1,
  'checks': 1024,
  'dist_threshold': 70000,
# sift or chog
  'descriptor': 'sift',
# euclidean, manhattan, minkowski, hik, hellinger, cs, kl
  'distance_type': 'euclidean',
# use >1 for weighted
  'num_neighbors': 1,
# highest, ransac, top_n, matchonce
  'vote_method': 'highest',
# custom configuration notation
  'confstring': '',
}

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
  nn = params['num_neighbors']
  nn = '' if nn == 1 else (',nn=%d' % nn)
  vote_method = '' if params['vote_method'] == 'highest' else ',%s' % params['vote_method']
  conf = ''
  if params['confstring']:
    conf = ',%s' % params['confstring']
  return '%s,threshold=%dk,searchparam=%d%s%s%s' % (indextype(params), params['dist_threshold']/1000, params['checks'], nn, vote_method, conf)

def run_parallel(dbdir, cells, querydir, querysift, outputFilePaths, params, num_threads=cpu_count()):
  semaphore = threading.Semaphore(num_threads)
  threads = []
  INFO("running queries with up to %d threads" % num_threads)
  for cell, outputFilePath in zip(cells, outputFilePaths):
    thread = Query(dbdir, cell, querydir, querysift, outputFilePath, params, semaphore)
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
  def __init__(self, celldir, cell, qdir, qfile, outfile, params=PARAMS_DEFAULT, barrier=None):
    threading.Thread.__init__(self)
    self.qpath = qdir + qfile
    self.cellpath = celldir + cell
    self.celldir = celldir
    self.outfile = outfile
    self.params = params
    self.barrier = barrier
    self.dump = self.outfile + '-detailed.npy'
    pyflann.set_distance_type(params['distance_type'])
    self.reader = reader.get_reader(params['descriptor'])

  def run(self):
    if self.barrier:
      self.barrier.acquire()
    if os.path.exists(self.outfile) and os.path.exists(self.dump):
      return
    self.flann = pyflann.FLANN()
    start = time.time()
    dataset, mapping = self._build_index()
    queryset = self.reader.load_file(self.qpath)
    qtime = time.time()
    results, dists = self.flann.nn_index(queryset['vec'], **self.params)
    INFO_TIMING("query took %f seconds" % (time.time() - qtime))
    sorted_counts = self.vote(queryset, dataset, mapping, results, dists)
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
      'top_n': self._vote_top_n,
      'ransac': self._vote_ransac,
    }[self.params['vote_method']](queryset, dataset, mapping, results, dists)
    return counts

  def _vote_top_n(self, queryset, dataset, mapping, results, dists):
    """Like vote_matchonce, but up to 3 db images per query feature.
       requires more than 1 nearest neighbor for results.
       Note that each image gets 1 vote max."""
    accept, reject, matchonce = 0, 0, 0
    counts = {} # map from img to counts
    closed = set()
    for i, dist_array in enumerate(dists):
      best = dist_array[0]
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
