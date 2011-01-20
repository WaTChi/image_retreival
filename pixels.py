# Find real locations of corresponding SIFT features.
# Provides map from (SIFTFile, pixel) => (lat, lon, alt)

from config import *
import numpy as np
import os

class LazyPixelMap:
  def __init__(self, filedir):
    self.filedir = filedir
    self.datastore = os.path.join(os.path.parent(filedir), 'pixeldata')
    if not os.path.isdir(datastore):
      os.mkdir(datastore)

  def ddFetch(siftfile, view):
    "Fetches pixels for view using ddObject"
    info = eval(open(os.path.join(self.filedir, view + '.info')).read())
    pixels = set()
    for x,y,foo,bar in SIFTReader.load_file(siftfile)['geom']:
      pixels.add((int(x), int(y)))
    data = ddGetAllPixels(pixels, info['id'], keep_None=True)
    assert len(data) == len(pixels)
    return data
  
  def open(siftfile):
    """Returns map of (x,y) => (lat, lon, alt)
       Assumes first 32 chars of filename are unique view"""
    view = os.path.basename(siftfile)[:32]
    cached = os.path.join(self.datastore, os.path.basename(view))
    if not os.path.exists(cached):
      data = self.ddFetch(siftfile, view)
      save_atomic(lambda d: np.save(d, data), cached)
    return np.load(data).item()
