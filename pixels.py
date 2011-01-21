# Find real locations of corresponding features.
# Provides map from (file, pixel) => (lat, lon, alt)
# Use pixelmap.open(...) for efficency.

from config import *
from reader import get_reader
from earthMine import ddGetAllPixels
import numpy as np
import os

class LazyPixelMap:
  def __init__(self, infodir):
    self.infodir = infodir
    self.datastore = os.path.join(os.path.dirname(infodir), 'pixeldata')
    if not os.path.isdir(self.datastore):
      os.mkdir(self.datastore)

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
  
  def open(self, featurefile):
    """Returns map of (x,y) => (lat, lon, alt)"""
    view = os.path.basename(featurefile)[:-8]
    cached = os.path.join(self.datastore, view) + '.npy'
    if not os.path.exists(cached):
      data = self.ddFetch(featurefile, view)
      save_atomic(lambda d: np.save(d, data), cached)
    return np.load(cached).item()
