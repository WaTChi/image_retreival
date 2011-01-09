#!/usr/bin/python
#
# Loads sift files, caching feature vectors for future runs.
# Use load_file(sift.txt) to read features from a single file.
# Use npy_cached_load(directory) to load a cell.
#

from config import *
import numpy as np
import pickle
import os
import shutil

NUM_DIMENSIONS = 128
IS_SIFT = lambda filename: 'sift.txt' in filename

# 'r' to enable numpy mmap loading
# it's marginally faster but occasionally hangs flann in 'D' state??
MMAP_MODE = None

def sift_iterator(siftname):
  """Returns feature values in chunks of arbitrary size."""
  with open(siftname) as data:
    count = -2
    for line in data:
      count += 1
      if count < 0 or count % 8 == 0:
        continue
      yield np.fromstring(line, sep=' ', dtype=np.uint8)

def write_features_to_ndarray(siftname, offset, dataset, key=None, keyset=None):
  """Adds features from a sift file to a dataset.
     Returns the new offset into the matrix."""
  step = 0
  for chunk in sift_iterator(siftname):
    # copy chunk into array
    dataset[offset][step:step+len(chunk)] = chunk
    step += len(chunk)
    if step >= NUM_DIMENSIONS:
      step = 0
      if key is not None:
        keyset[offset] = key
      offset += 1
  return offset

def load_file(siftname):
  """Loads single sift file into numpy array."""
  with open(siftname) as f:
    num_features = int(f.readline().split()[0])
  dataset = np.ndarray((num_features, NUM_DIMENSIONS), dtype=np.uint8)
  write_features_to_ndarray(siftname, 0, dataset)
  return dataset

def load_file_lossy(siftname, ratio):
  """Simulates quantization compression of SIFT vectors. Ratio should be power of 2"""
  dataset = load_file(siftname)
  assert ratio == int(ratio) and ratio >= 1
  if ratio != 1:
    unit = 256/2**(8/ratio)
    dataset = dataset/unit*unit + unit/2
  return dataset

def npy_save_sift_directory(directory, cellid):
  """Writes all sift features found in a directory to a file.
     Also builds a reverse lookup table."""
  num_features = 0
  # first count number of total features
  for name in os.listdir(directory):
    if IS_SIFT(name):
      with open(os.path.join(directory, name)) as f:
        num_features += int(f.readline().split()[0])
  dataset = np.ndarray((num_features, NUM_DIMENSIONS), dtype=np.uint8)
  keyset = np.ndarray((num_features, 1), dtype=np.uint16)
  offset = 0
  key = 0
  lookup_table = {}
  # now begin the actual read
  for name in os.listdir(directory):
    if IS_SIFT(name):
      offset = write_features_to_ndarray(os.path.join(directory, name), offset, dataset, key, keyset)
      lookup_table[key] = name
      if key % 200 == 0:
        INFO('%d/%d features read' % (offset, num_features))
      key += 1
  for dest in getdests(directory, cellid + '-features.npy'):
    save_atomic(lambda d: np.save(d, dataset), dest)
  for dest in getdests(directory, cellid + '-keys.npy'):
    save_atomic(lambda d: np.save(d, keyset), dest)
  for dest in getdests(directory, cellid + '-map.p'):
    save_atomic(lambda d: pickle.dump(lookup_table, open(d, 'w')), dest)

def getdests(directory, name):
  default = os.path.join(os.path.dirname(directory), name)
  if IS_REMOTE(directory):
    INFO('copying %s to local cache' % name)
    local = os.path.join(CACHE_PATH, name)
    return [default, local]
  return [default]

def getfile(directory, name):
  default = os.path.join(os.path.dirname(directory), name)
  if IS_REMOTE(directory):
    local = os.path.join(CACHE_PATH, name)
    if os.path.exists(local):
      return local
    elif os.path.exists(default):
      INFO('copying %s to local cache' % name)
      save_atomic(lambda d: shutil.copyfile(default, d), local)
      return local
  return default

def getcellid(directory):
  return os.path.basename(directory)

def npy_cached_load(directory):
  """Efficiently loads a matrix of sift features and reverse lookup table
     for a directory of sift files."""
  cellid = getcellid(directory)
  data_out = getfile(directory, cellid + '-features.npy')
  key_out = getfile(directory, cellid + '-keys.npy')
  map_out = getfile(directory, cellid + '-map.p')
  if not os.path.exists(data_out) or not os.path.exists(map_out) or not os.path.exists(key_out):
    INFO('finding sift files')
    npy_save_sift_directory(directory, cellid)
  mapping = pickle.load(open(map_out))
  return np.load(data_out, mmap_mode=MMAP_MODE), mapping, np.load(key_out, mmap_mode=MMAP_MODE)

if __name__ == '__main__':
  from sys import argv
  data, mapping, keyset = npy_cached_load(argv[1])
  INFO('loaded %d features from %d images' % (len(data), len(mapping)))

# vim: et sw=2
