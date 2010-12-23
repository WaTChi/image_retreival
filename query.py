#!/usr/bin/python

from config import *
from SIFTReader import *
import pyflann
import numpy as np
import os

PARAMS_DEFAULT = {
  'algorithm': 'kdtree1',
  'checks': 128,
  'nn': 5,
}

class Query:
  def __init__(self, celldir, cell, qdir, qfile, outfile, params=PARAMS_DEFAULT):
    self.qpath = qdir + qfile
    self.cellpath = celldir + cell
    self.outfile = outfile
    self.params = params

  def run(self):
    records, mapping = npy_cached_load(self.cellpath)
    query = load_file(self.qpath)
    queryset = FEATURE_VIEW(query)
    flann = pyflann.FLANN()
    params = flann.build_index(FEATURE_VIEW(records), kwargs=self.params)
    INFO(params)
    results, dists = flann.nn_index(queryset, self.params['nn'], checks=self.params['checks'])
    INFO(results)
    INFO(dists)

# vim: et sw=2
