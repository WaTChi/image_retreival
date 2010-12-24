#!/usr/bin/python

from config import *
from SIFTReader import *
import pyflann
import numpy as np
import os

PARAMS_DEFAULT = {
  'algorithm': 'kdtree',
  'trees': 1,
  'checks': 128,
  'num_neighbors': 1,
  'dist_threshold': 70000,
}

def str_params(params):
  return ','.join(sorted(['%s=%s' % (k,v) for k,v in params.iteritems()]))

class Query:
  def __init__(self, celldir, cell, qdir, qfile, outfile, params=PARAMS_DEFAULT):
    self.qpath = qdir + qfile
    self.cellpath = celldir + cell
    self.outfile = outfile
    self.params = params
    self.flann = pyflann.FLANN()

  def _build_index(self):
    iname = getcellid(self.cellpath) + '-' + self.params['algorithm'] + str(self.params.get('trees', '')) + '.uint8.index'
    index = getfile(self.cellpath, iname)
    dataset, mapping, keyset = npy_cached_load(self.cellpath)
    if os.path.exists(index):
      INFO('loading index %s' % iname)
      self.flann.load_index(index, dataset)
      return mapping, keyset
    INFO('creating %s' % iname)
    INFO(self.flann.build_index(dataset, algorithm=self.params['algorithm'], trees=self.params.get('trees', 0)))
    for out in getdests(self.cellpath, iname):
      self.flann.save_index(out)
    return mapping, keyset

  def run(self):
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
    with open(self.outfile + '.test', 'w') as f:
      for tally in votes:
        f.write("%d\t%s\n" % tally)
        total += tally[0]
    INFO('put %d/%d votes into %d bins' % (total, len(results), len(votes)))

# vim: et sw=2
