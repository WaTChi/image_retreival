#!/usr/bin/python

from config import *
from SIFTReader import *
import pyflann
import os

PARAMS_DEFAULT = {
  'algorithm': 'kdtree',
  'trees': 1,
  'checks': 1024,
  'num_neighbors': 1,
  'log_level': 'info',
# for convenience; not a flann build parameter
  'dist_threshold': 70000,
  'distance_type': 'euclidean',
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
    INFO("running query")
    queryset = load_file(self.qpath)
    counts = {} # map from sift index to counts
    results, dists = self.flann.nn_index(queryset, **self.params)
    for i, dist in enumerate(dists):
      if dist < self.params['dist_threshold']:
        k = keyset[results[i]][0]
        if k in counts:
          counts[k] += 1
        else:
          counts[k] = 1
    votes = sorted([(counts[index], mapping[index]) for index in counts], reverse=True)
    total = 0
    with open(self.outfile, 'w') as f:
      for tally in votes:
        f.write("%d\t%s\n" % tally)
        total += tally[0]
    INFO("took %f total" % (time.time() - start))
    INFO('put %d/%d votes into %d bins' % (total, len(results), len(votes)))

  def _build_index(self):
    iname = '%s-%s.uint8.index' % (getcellid(self.cellpath), indextype(self.params))
    index = getfile(self.cellpath, iname)
    start = time.time()
    dataset, mapping, keyset = npy_cached_load(self.cellpath)
    INFO("load took %f seconds" % (start - time.time()) )
    if os.path.exists(index):
      INFO('loading index %s' % iname)
      self.flann.load_index(index, dataset)
      return mapping, keyset
    INFO('creating %s' % iname)
    INFO(self.flann.build_index(dataset, **self.params))
    for out in getdests(self.cellpath, iname):
      self.flann.save_index(out)
    return mapping, keyset

# vim: et sw=2
