#!/usr/bin/env python

# SURF extractor script

import numpy as np
import os
import time
import cv
from reader import SURFReader
from config import *

def extract_surf(jpgfile):
  start = time.time()
  out = os.path.join(os.path.dirname(jpgfile), os.path.basename(jpgfile)[:-4] + 'surf.npy')
  if os.path.exists(out):
    INFO('%s already exists' % out)
    return

  im = cv.LoadImageM(jpgfile, cv.CV_LOAD_IMAGE_GRAYSCALE)
  INFO('cv loaded %dx%d image' % (im.rows, im.cols))

  g, features = cv.ExtractSURF(im, None, cv.CreateMemStorage(), (0, 500, 3, 4))
  data = np.ndarray(len(features), SURFReader.surf_dtype)

  for i in range(len(features)):
    data[i]['vec'] = np.fromiter(features[i], np.float32)
    data[i]['geom'] = np.fromiter([g[i][0][0], g[i][0][1], g[i][2]], np.uint16)
    data[i]['index'] = 0

## Simple Quantization into bytes
#  for i in range(len(features)):
#    surfvalues = np.fromiter(features[i], np.float)
#
#    assert max(surfvalues) <= 1.0
#    assert min(surfvalues) >= -1.0
#
#    data[i]['vec'] = np.int8(127*surfvalues)
#    data[i]['geom'] = np.fromiter([g[i][0][0], g[i][0][1], g[i][2]], np.uint16)
#    data[i]['index'] = 0

  save_atomic(lambda d: np.save(d, data), out)
  INFO('cv wrote %d features' % len(features))
  INFO_TIMING('took %f seconds' % (time.time() - start))

if __name__ == '__main__':
  import sys
  INFO('extracting SURF features from jpgs in %s' % sys.argv[1])
  for file in os.listdir(sys.argv[1]):
    if file.lower().endswith('.jpg'):
      extract_surf(os.path.join(sys.argv[1], file))

# vim: et sw=2
