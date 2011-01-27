#!/usr/bin/python
#
# Finds nearest neighbors of a set of features in a database.
# Use Query(...).run() for a single query
# or
# run_parallel(...) to query multiple cells in parallel
#

from config import *
import reader
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
# highest, weighted, revote_exact, location
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
  vote_method = '' if params['vote_method'] == 'highest' else ',%s' % params['vote_method']
  conf = ''
  if params['confstring']:
    conf = ',%s' % params['confstring']
  return '%s,threshold=%dk,searchparam=%d%s%s' % (indextype(params), params['dist_threshold']/1000, params['checks'], vote_method, conf)

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
  def __init__(self, celldir, cell, qdir, qfile, outfile, params=PARAMS_DEFAULT, barrier=None, dump=None):
# for more detailed results, specify dump.
# dump will be a *.npy file with:
# map of (imagename => list of (vector from db, matching vector from query))
    threading.Thread.__init__(self)
    self.qpath = qdir + qfile
    self.cellpath = celldir + cell
    self.outfile = outfile
    self.params = params
    self.barrier = barrier
    self.dump = dump
    pyflann.set_distance_type(params['distance_type'])
    self.reader = reader.get_reader(params['descriptor'])

  def run(self):
    if self.barrier:
      self.barrier.acquire()
    self.flann = pyflann.FLANN()
    start = time.time()
    dataset, mapping = self._build_index()
    queryset = self.reader.load_file(self.qpath)
    qtime = time.time()
    results, dists = self.flann.nn_index(queryset['vec'], **self.params)
    INFO_TIMING("query took %f seconds" % (time.time() - qtime))
    counts = self.vote(queryset, dataset, mapping, results, dists)
    if self.dump:
      INFO('saving raw vote data in %s' % self.dump)
      np.save(self.dump, counts)
    votes = sorted([(len(counts[img]), img) for img in counts], reverse=True)
#    total = 0
    with open(self.outfile, 'w') as f:
      for tally in votes:
        f.write("%f\t%s\n" % tally)
#        total += tally[0]
    INFO_TIMING("took %f total" % (time.time() - start))
    if self.barrier:
      self.barrier.release()
    # release memory - in case the Query is still around
    self.flann = None

  def vote(self, queryset, dataset, mapping, results, dists):
    INFO('voting with method %s' % self.params['vote_method'])
    counts = {
      'highest': self._vote_highest,
      'weighted': self._vote_weighted,
      'location': self._vote_location,
      'revote_exact': self._vote_revote_exact,
    }[self.params['vote_method']](queryset, dataset, mapping, results, dists)
    return counts

  def _vote_highest(self, queryset, dataset, mapping, results, dists):
    counts = {} # map from sift index to counts
    accept, reject = 0, 0
    for i, dist in enumerate(dists):
      if dist < self.params['dist_threshold']:
        accept += 1
        img = mapping[dataset[results[i]]['index']]
        if img not in counts:
          counts[img] = []
        counts[img].append(results[i])
      else:
        reject += 1
    INFO('accepted %d/%d votes' % (accept, accept + reject))
    return counts

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
