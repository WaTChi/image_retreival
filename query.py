#!/usr/bin/python
#
# Finds nearest neighbors of a set of features in a database.
# Use Query(...).run() for a single query
# or
# run_parallel(...) to query multiple cells in parallel
#
# AUTHOR: Eric Liang

from config import *
from SIFTReader import *
from multiprocessing import cpu_count
import info
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
# euclidean, manhattan, minkowski, hik, hellinger, cs, kl
  'distance_type': 'euclidean',
# use >1 for weighted
  'num_neighbors': 1,
# highest, weighted, revote_exact, location
  'vote_method': 'highest',
}

# I'm not actually sure if the distance function affects
# index compatibility. Someone check please?
def indextype(params):
  alg = params['algorithm']
  dtype = params['distance_type']
  distname = '' if dtype == 'euclidean' else ('-%s' % dtype)
  if alg == 'kdtree':
    return 'kdtree%d%s' % (params['trees'], distname)
  return '%s%s' % (alg, distname)

def searchtype(params):
  vote_method = '' if params['vote_method'] == 'highest' else ',%s' % params['vote_method']
  return '%s,threshold=%dk,searchparam=%d%s' % (indextype(params), params['dist_threshold']/1000, params['checks'], vote_method)

def run_parallel(dbdir, cells, querydir, querysift, outputFilePaths, params, num_threads=cpu_count()):
  semaphore = threading.Semaphore(num_threads)
  threads = []
  ok = False
  for cell, outputFilePath in zip(cells, outputFilePaths):
    if not os.path.exists(outputFilePath):
      if not ok:
        INFO("running queries with %d threads" % num_threads)
        ok = True
      thread = Query(dbdir, cell, querydir, querysift, outputFilePath, params, semaphore)
      threads.append(thread)
      thread.start()
  for thread in threads:
     thread.join()

class Query(threading.Thread):
  def __init__(self, celldir, cell, qdir, qfile, outfile, params=PARAMS_DEFAULT, barrier=None):
    threading.Thread.__init__(self)
    self.qpath = qdir + qfile
    self.cellpath = celldir + cell
    self.outfile = outfile
    self.params = params
    self.barrier = barrier
    pyflann.set_distance_type(params['distance_type'])
    self.flann = pyflann.FLANN()

  def run(self):
    if self.barrier:
      self.barrier.acquire()
    start = time.time()
    mapping, keyset = self._build_index()
    queryset = load_file(self.qpath)
    qtime = time.time()
    results, dists = self.flann.nn_index(queryset, **self.params)
    INFO_TIMING("query took %f seconds" % (time.time() - qtime))
    votes = self.vote(queryset, mapping, keyset, results, dists)
    total = 0
    with open(self.outfile, 'w') as f:
      for tally in votes:
        f.write("%f\t%s\n" % tally)
        total += tally[0]
    INFO_TIMING("took %f total" % (time.time() - start))
    if self.barrier:
      self.barrier.release()

  def vote(self, queryset, mapping, keyset, results, dists):
    INFO('voting with method %s' % self.params['vote_method'])
    return {
      'highest': self._vote_highest,
      'weighted': self._vote_weighted,
      'location': self._vote_location,
      'revote_exact': self._vote_revote_exact,
    }[self.params['vote_method']](queryset, mapping, keyset, results, dists)

  def _vote_highest(self, queryset, mapping, keyset, results, dists):
    counts = {} # map from sift index to counts
    for i, dist in enumerate(dists):
      if dist < self.params['dist_threshold']:
        k = keyset[results[i]][0]
        counts[k] = 1.0 + counts.get(k, 0.0)
    votes = sorted([(counts[index], mapping[index]) for index in counts], reverse=True)
    return votes

  def extract(self, dataset, keyset, indices):
    "returns set of features from images in indices"
    indices = set(indices)
    condition = np.zeros(len(keyset), np.bool)
    for i in indices:
      condition |= (keyset.reshape(-1) == i)
    sdataset = np.compress(condition, dataset, axis=0)
    skeyset = np.compress(condition, keyset, axis=0)
    return sdataset, skeyset

  def _vote_revote_exact(self, queryset, mapping, keyset, results, dists):
    counts = {} # map from sift index to counts
    for i, dist in enumerate(dists):
      if dist < self.params['dist_threshold']:
        k = keyset[results[i]][0]
        counts[k] = 1.0 + counts.get(k, 0.0)
    # list of (votes, index_of_image)
    top_ten = sorted([(counts[index], index) for index in counts], reverse=True)[:5]
    sdataset, skeyset = self.extract(self.dataset, keyset, [i for (v, i) in top_ten])
    results, dists = self.flann.nn(sdataset, queryset, 1, algorithm='linear')
    votes = self._vote_highest(queryset, mapping, skeyset, results, dists)
    return votes

  def _vote_weighted(self, queryset, mapping, keyset, results, dists):
    """Keep voting for nearest N, weighting by (dist/MIN_DISTANCE)^-k
       requires more than 1 nearest neighbor for results.
       Note that each image gets 1 vote max."""
    k = 10.0
    accept, reject = 0, 0
    counts = {} # map from sift index to counts
    for i, dist_array in enumerate(dists):
      best = dist_array[0]
      marked = set()
      for j, dist in enumerate(dist_array):
        if dist < self.params['dist_threshold']:
          image = keyset[results[i][j]][0]
          if image not in marked:
            counts[image] = 1.0*(dist/best)**-k + counts.get(image, 0.0)
            marked.add(image)
          accept += 1
        else:
          reject += 1
    INFO('threw away %d/%d votes' % (reject, accept + reject))
    votes = sorted([(counts[index], mapping[index]) for index in counts], reverse=True)
    return votes

  def _vote_location(self, queryset, mapping, keyset, results, dists):
    TOP_N = 10 # O(N^2)
    counts = {} # map from sift index to counts
    for i, dist in enumerate(dists):
      if dist < self.params['dist_threshold']:
        k = keyset[results[i]][0]
        counts[k] = 1.0 + counts.get(k, 0.0)
    # get top N of counts
    counts = dict([(k,v) for (v,k) in sorted([(v,k) for (k,v) in counts.iteritems()], reverse=True)[:TOP_N]])
    w = {} # weight by exponential dropoff of distance
    for source_index in counts:
      for test_index in counts:
        dist = info.siftdistance(mapping[source_index], mapping[test_index])
        # 1.10^-dist gives 40% at 10 meters, 15% at 20 meters
        # 100 is infinitely steep dropoff
        w[test_index] = w.get(test_index, 0.0) + counts[source_index]*(1.10**(-dist))
    votes = sorted([(w[index], mapping[index]) for index in w], reverse=True)
    return votes

  def _build_index(self):
    start = time.time()
    iname = '%s-%s.uint8.index' % (getcellid(self.cellpath), indextype(self.params))
    index = getfile(self.cellpath, iname)
    dataset, mapping, keyset = npy_cached_load(self.cellpath)
    self.dataset = dataset
    INFO_TIMING("dataset load took %f seconds" % (time.time() - start))
    if os.path.exists(index):
      s = time.time()
      self.flann.load_index(index, dataset)
      INFO_TIMING("index load took %f seconds" % (time.time() - s))
      return mapping, keyset
    INFO('creating %s' % iname)
    start = time.time()
    INFO(self.flann.build_index(dataset, **self.params))
    INFO_TIMING("index creation took %f seconds" % (time.time() - start))
    for out in getdests(self.cellpath, iname):
      self.flann.save_index(out + '.tmp')
      os.rename(out + '.tmp', out)
    return mapping, keyset

# vim: et sw=2
