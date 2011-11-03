#!/usr/bin/env python
# Find real locations of corresponding features.
# Use pixelmap.open(siftfile) for reading a specific file.
# Use pixelmap.get(siftfile, x, y) for general fast reads.
# Provides map from (pixel) => (lat, lon, alt)

import time
from config import *
from reader import get_reader
from earthMine import ddGetAllPixels
import numpy as np
import os

class PixelMap:
  def __init__(self, infodir):
    self.infodir = infodir
    self.datastore = os.path.join(os.path.dirname(os.path.abspath(infodir)), 'pixeldata')
    if not os.path.isdir(self.datastore):
      os.mkdir(self.datastore)
    self.cached = {}

  def ddFetch(self, featurefile, view):
    "Fetches pixels for view using ddObject"
    info = eval(open(os.path.join(self.infodir, view + '.info')).read())
    pixels = set()
    reader = get_reader(os.path.basename(featurefile))
    for x,y,foo,bar in reader.load_file(featurefile)['geom']:
      pixels.add((int(x), int(y)))
    data = ddGetAllPixels(pixels, info['id'], keep_None=True)
    assert len(data) == len(pixels)
    return data
  
  def get(self, featurefile, x, y):
    start = time.time()
    x, y = int(x), int(y)
    if featurefile in self.cached:
      pixmap = self.cached[featurefile]
    else:
      pixmap = self.cached[featurefile] = self.open(featurefile)
    return pixmap[x,y]
  
  def open(self, featurefile):
    """Returns map of (x,y) => {'lat':lat, lon, alt}"""
    name = os.path.basename(featurefile)[:-4] # gps coords, angle
    view = os.path.basename(featurefile)[:-8] # + descriptor type
    cached = os.path.join(self.datastore, name) + '.npy'
    if not os.path.exists(cached):
      INFO("*** fetching pixel data from earthmine ***")
      data = self.ddFetch(featurefile, view)
      save_atomic(lambda d: np.save(d, data), cached)
    return np.load(cached).item()

if __name__ == '__main__':
  mapper = PixelMap('/home/ericl/shiraz/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829')
  superdir = '/home/ericl/shiraz/Research/cells/g=100,r=d=236.6/'
  for dir in os.listdir(superdir):
    dir = os.path.join(superdir, dir)
    if os.path.isdir(dir):
      for f in get_reader('sift').get_feature_files_in_dir(dir):
        try:
          mapper.open(f)
        except:
          pass

# vim: et sw=2
