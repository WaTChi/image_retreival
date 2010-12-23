#!/usr/bin/python

from config import *
import numpy as np
import pickle
import os
import shutil

NUM_DIMENSIONS = 128
ROW = np.dtype([('image', np.uint16), ('features', np.uint8, (NUM_DIMENSIONS,))])
FEATURE_VIEW = lambda dataset: dataset.getfield(ROW['features'], 2)
IMAGE_VIEW = lambda dataset: dataset.getfield(ROW['image'], 0)
IS_SIFT = lambda filename: 'sift.txt' in filename

def INFO(x):
  print 'INFO: ' + str(x)

def sift_iterator(siftname):
  """Returns feature values in chunks of arbitrary size."""
  with open(siftname) as data:
    count = -2
    for line in data:
      count += 1
      if count < 0 or count % 8 == 0:
        continue
      yield np.fromstring(line, sep=' ', dtype=np.uint8)

def write_features_to_ndarray(siftname, offset, dataset, key):
  """Adds features from a sift file to a dataset.
     Returns the new offset into the matrix."""
  dataset[offset][0][0] = key
  step = 0
  for chunk in sift_iterator(siftname):
    # copy chunk into array
    dataset[offset][0][1][step:step+len(chunk)] = chunk
    step += len(chunk)
    if step >= NUM_DIMENSIONS:
      step = 0
      offset += 1
  return offset

def npy_save_sift_directory(directory, outname):
  """Writes all sift features found in a directory to a file.
     Also builds a reverse lookup table."""
  num_features = 0
  # first count number of total features
  for name in os.listdir(directory):
    if IS_SIFT(name):
      with open(os.path.join(directory, name)) as f:
        num_features += int(f.readline().split()[0])
  lookup_table = {}
  dataset = np.ndarray((num_features, 1), dtype=ROW)
  offset = 0
  key = 0
  lookup_table = {}
  # now begin the actual read
  for name in os.listdir(directory):
    if IS_SIFT(name):
      offset = write_features_to_ndarray(os.path.join(directory, name), offset, dataset, key)
      lookup_table[key] = name
      key += 1
      INFO('%d/%d features read' % (offset, num_features))
  for dest in getdests(directory, outname + '.npy'):
    np.save(dest, dataset)
  for dest in getdests(directory, outname + '.map'):
    pickle.dump(lookup_table, open(dest, 'w'))

def getdests(directory, name):
  default = os.path.join(directory, name)
  if IS_REMOTE(directory):
    INFO('copying %s to local cache' % name)
    local = os.path.join(CACHE_PATH, name)
    return [default, local]
  return [default]

def findfile(directory, name):
  default = os.path.join(directory, name)
  if IS_REMOTE(directory):
    local = os.path.join(CACHE_PATH, name)
    if os.path.exists(local):
      INFO('using local copy of %s' % name)
      return local
    elif os.path.exists(default):
      INFO('copying %s to local cache' % name)
      shutil.copyfile(default, local)
      return local
  return default

def npy_cached_load(directory, map_only=False):
  """Efficiently loads a matrix of sift features and reverse lookup table
     for a directory of sift files."""
  names = []
  for name in os.listdir(directory):
    if IS_SIFT(name):
      names.append(name)
  names.sort()
  outname = 'data-' + hex(hash(tuple(names)))[2:]
  data_out = findfile(directory, outname + '.npy')
  map_out = findfile(directory, outname + '.map')
  if not os.path.exists(data_out) or not os.path.exists(map_out):
    INFO('finding sift files')
    npy_save_sift_directory(directory, outname)
    data_out = findfile(directory, outname + '.npy')
    map_out = findfile(directory, outname + '.map')
  record = {'mapping': pickle.load(open(map_out))}
  if not map_only:
    record['data'] = np.load(data_out)
  return record

if __name__ == '__main__':
  from sys import argv
  data = npy_cached_load(argv[1])
  INFO('loaded %d features from %d images' % (len(data['data']), len(data['mapping'])))

# vim: et sw=2
