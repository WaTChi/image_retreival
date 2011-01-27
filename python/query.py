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
    INFO('saving detailed vote data in %s' % self.dump)
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
      'ransac': self._vote_ransac,
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
        counts[img].append({'db': dataset[results[i]]['geom'].copy(),
                            'query': queryset[i]['geom'].copy()})
      else:
        reject += 1
    INFO('accepted %d/%d votes' % (accept, accept + reject))
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    return sorted_counts

  def _vote_ransac(self, queryset, dataset, mapping, results, dists):
    sorted_counts = self._vote_highest(queryset, dataset, mapping, results, dists)
    # filters out outliers from counts until
    # filtered_votes(ith image) > votes(jth image) for all j != i
    # and returns top 10 filtered results
    raise NotImplementedError

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
