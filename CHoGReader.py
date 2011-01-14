#!/usr/bin/python
#
# Loads chog files, caching feature vectors for future runs.
# Use load_file(chog.txt) to read features from a single file.
# Use npy_cached_load(directory) to load a cell.
#

from config import *
import numpy as np
import pickle
import os
import shutil

NUM_DIMENSIONS = 63
IS_CHOG = lambda filename: 'chog.txt' in filename

def chog_iterator(chogname):
  """Returns feature values in chunks of NUM_DIMENSIONS."""
  with open(chogname) as data:
    for line in data:
      yield np.fromstring(line, sep=' ', dtype=np.float)[5:]

def write_features_to_ndarray(chogname, offset, dataset, key=None, keyset=None):
  """Adds features from a chog file to a dataset.
     Returns the new offset into the matrix."""
  for chunk in chog_iterator(chogname):
    # copy chunk into array
    assert len(chunk) == NUM_DIMENSIONS
    dataset[offset] = chunk
    if key is not None:
      keyset[offset] = key
    offset += 1
  return offset

def load_file(chogname):
  """Loads single chog file into numpy array."""
  with open(chogname) as f:
    num_features = f.read().count('\n')
  dataset = np.ndarray((num_features, NUM_DIMENSIONS), dtype=np.uint8)
  write_features_to_ndarray(chogname, 0, dataset)
  return dataset

def npy_save_chog_directory(directory, cellid):
  """Writes all chog features found in a directory to a file.
     Also builds a reverse lookup table."""
  num_features = 0
  # first count number of total features
  for name in os.listdir(directory):
    if IS_CHOG(name):
      with open(os.path.join(directory, name)) as f:
        num_features += int(f.readline().split()[0])
  dataset = np.ndarray((num_features, NUM_DIMENSIONS), dtype=np.uint8)
  keyset = np.ndarray((num_features, 1), dtype=np.uint16)
  offset = 0
  key = 0
  lookup_table = {}
  # now begin the actual read
  for name in os.listdir(directory):
    if IS_CHOG(name):
      offset = write_features_to_ndarray(os.path.join(directory, name), offset, dataset, key, keyset)
      lookup_table[key] = name
      if key % 200 == 0:
        INFO('%d/%d features read' % (offset, num_features))
      key += 1
  for dest in getdests(directory, cellid + '-features-chog.npy'):
    save_atomic(lambda d: np.save(d, dataset), dest)
  for dest in getdests(directory, cellid + '-keys-chog.npy'):
    save_atomic(lambda d: np.save(d, keyset), dest)
  for dest in getdests(directory, cellid + '-map-chog.p'):
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
  """Efficiently loads a matrix of chog features and reverse lookup table
     for a directory of chog files."""
  cellid = getcellid(directory)
  data_out = getfile(directory, cellid + '-features-chog.npy')
  key_out = getfile(directory, cellid + '-keys-chog.npy')
  map_out = getfile(directory, cellid + '-map-chog.p')
  if not os.path.exists(data_out) or not os.path.exists(map_out) or not os.path.exists(key_out):
    INFO('finding chog files')
    npy_save_chog_directory(directory, cellid)
  mapping = pickle.load(open(map_out))
  return np.load(data_out), mapping, np.load(key_out)

if __name__ == '__main__':
  from sys import argv
  data, mapping, keyset = npy_cached_load(argv[1])
  INFO('loaded %d features from %d images' % (len(data), len(mapping)))

# vim: et sw=2
