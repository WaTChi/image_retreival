#!/usr/bin/python

from config import *
from SIFTReader import *
import info
import time
import pyflann
import os

PARAMS_DEFAULT = {
  'algorithm': 'kdtree',
  'trees': 1,
  'checks': 1024,
# for convenience; not flann build parameters
  'dist_threshold': 70000,
  'distance_type': 'euclidean',
# highest, threshold, exp_dropoff
  'vote_method': 'highest',
# set nn >1 for threshold vote method
  'num_neighbors': 1,
}

def siftdistance(a, b):
  if a == b:
    return 0
  lat1, lon1 = info.getSIFTCoord(a)
  lat2, lon2 = info.getSIFTCoord(b)
  return info.distance(lat1, lon1, lat2, lon2)

# I'm not actually sure if the distance function affects
# index compatibility. Someone check please?
def indextype(params):
  alg = params['algorithm']
  dtype = params['distance_type']
  distname = '' if dtype == 'euclidean' else ('-%s' % dtype)
  if alg == 'kdtree':
    return 'kdtree%d%s' % (params['trees'], distname)
  return '%s%s' % (alg, distname)

class Query:
  def __init__(self, celldir, cell, qdir, qfile, outfile, params=PARAMS_DEFAULT):
    self.qpath = qdir + qfile
    self.cellpath = celldir + cell
    self.outfile = outfile
    self.params = params
    pyflann.set_distance_type(params['distance_type'])
    self.flann = pyflann.FLANN()

  def run(self):
    start = time.time()
    mapping, keyset = self._build_index()
    queryset = load_file(self.qpath)
    qtime = time.time()
    results, dists = self.flann.nn_index(queryset, **self.params)
    INFO("query took %f seconds" % (time.time() - qtime))
    votes = self.vote(queryset, mapping, keyset, results, dists)
    total = 0
    with open(self.outfile, 'w') as f:
      for tally in votes:
        f.write("%d\t%s\n" % tally)
        total += tally[0]
    INFO("took %f total" % (time.time() - start))

  def vote(self, queryset, mapping, keyset, results, dists):
    INFO('voting with method %s' % self.params['vote_method'])
    return {
      'highest': self._vote_highest,
      'exp_dropoff': self._vote_exp_dropoff,
      'threshold': self._vote_threshold,
    }[self.params['vote_method']](queryset, mapping, keyset, results, dists)

  def _vote_highest(self, queryset, mapping, keyset, results, dists):
    counts = {} # map from sift index to counts
    for i, dist in enumerate(dists):
      if dist < self.params['dist_threshold']:
        k = keyset[results[i]][0]
        if k in counts:
          counts[k] += 1
        else:
          counts[k] = 1
    votes = sorted([(counts[index], mapping[index]) for index in counts], reverse=True)
    return votes

  def _vote_exp_dropoff(self, queryset, mapping, keyset, results, dists):
    "idea is that we vote based on location, not the images"
    counts = {} # map from sift index to counts
    for i, dist in enumerate(dists):
      if dist < self.params['dist_threshold']:
        k = keyset[results[i]][0]
        if k in counts:
          counts[k] = 1.0 + counts.get(k, 0.0)
    w = {} # weight by exponential dropoff of distance
    for source_index in counts:
      for test_index in counts:
        dist = siftdistance(mapping[source_index], mapping[test_index])
        # 1.05^-dist gives 60% at 10 meters, 30% at 20 meters
        w[test_index] = w.get(test_index, 0.0) + counts[source_index]*(1.05^-dist)
    votes = sorted([(w[index], mapping[index]) for index in w], reverse=True)
    return votes

  def _vote_threshold(self, queryset, mapping, keyset, results, dists):
    "idea is that we allow multiple votes per feature if good enough"
    counts = {} # map from sift index to counts
    deviation_allowed = 3000
    for i, dist_array in enumerate(dists):
      best = self.params['dist_threshold'] - deviation_allowed
      for j, dist in enumerate(dist_array):
        if dist < best + deviation_allowed:
          best = dist
          k = keyset[results[i][j]][0]
          if k in counts:
            counts[k] = 1.0 + counts.get(k, 0.0)
    votes = sorted([(w[index], mapping[index]) for index in w], reverse=True)
    return votes

  def _build_index(self):
    start = time.time()
    iname = '%s-%s.uint8.index' % (getcellid(self.cellpath), indextype(self.params))
    index = getfile(self.cellpath, iname)
    dataset, mapping, keyset = npy_cached_load(self.cellpath)
    INFO("dataset load took %f seconds" % (time.time() - start))
    if os.path.exists(index):
      s = time.time()
      self.flann.load_index(index, dataset)
      INFO("index load took %f seconds" % (time.time() - s))
      return mapping, keyset
    INFO('creating %s' % iname)
    start = time.time()
    INFO(self.flann.build_index(dataset, **self.params))
    INFO("index creation took %f seconds" % (time.time() - start))
    for out in getdests(self.cellpath, iname):
      self.flann.save_index(out)
    return mapping, keyset

# vim: et sw=2
