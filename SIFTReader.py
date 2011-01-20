#!/usr/bin/python
#
# Loads sift files, caching feature vectors for future runs.
# Use load_file(sift.txt) to read features from a single file.
# Use npy_cached_load(directory) to load a cell.

from config import *
import numpy as np
import os

MMAP_MODE = 'r' # set to None if things hang
NUM_DIMENSIONS = 128
IS_SIFT = lambda filename: 'sift.txt' in filename

# I would use np.object instead of index but it costs ~500bits
DTYPE_ROW = {
  'names': ['vec', 'geom', 'index'],
  'formats': ['%duint8' % NUM_DIMENSIONS, '4float32', 'uint16'],
}

def sift_iterator(siftname):
  """Returns feature values in chunks of arbitrary size."""
  count = -2 # skip first line
  with open(siftname) as data:
    for line in data:
      count += 1
      if count < 0:
        continue
      if count % 8 == 0:
        yield np.fromstring(line, sep=' ', dtype=np.float32)
      else:
        yield np.fromstring(line, sep=' ', dtype=np.uint8)

def write_features_to_ndarray(siftname, offset, dataset, index):
  """Adds features from a sift file to a dataset.
     Returns the new offset into the matrix."""
  step = 0
  for chunk in sift_iterator(siftname):
    if chunk.dtype == np.float32:
      dataset['geom'][offset] = chunk
    else: # vector info
      dataset['vec'][offset][step:step+len(chunk)] = chunk
      step += len(chunk)
      if step >= NUM_DIMENSIONS:
        dataset['index'][offset] = index
        step = 0 # on to next vector
        offset += 1
  return offset

def load_file(siftname):
  """Loads single sift file into numpy array."""
  with open(siftname) as f:
    num_features = int(f.readline().split()[0])
    assert num_features > 0
  dataset = np.ndarray(num_features, DTYPE_ROW)
  write_features_to_ndarray(siftname, 0, dataset, 0)
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
  assert num_features > 0
  dataset = np.ndarray(num_features, DTYPE_ROW)
  offset = 0
  image_index = 0
  lookup_table = {}
  # now begin the actual read
  for name in os.listdir(directory):
    if IS_SIFT(name):
      offset = write_features_to_ndarray(os.path.join(directory, name), offset, dataset, image_index)
      lookup_table[image_index] = name
      if image_index % 200 == 0:
        INFO('%d/%d features read' % (offset, num_features))
      image_index += 1
  for dest in getdests(directory, cellid + '-sift.npy'):
    save_atomic(lambda d: np.save(d, dataset), dest)
  for dest in getdests(directory, cellid + '-pydata.npy'):
    save_atomic(lambda d: np.save(d, lookup_table), dest)

def npy_cached_load(directory):
  """Efficiently loads a matrix of sift features and reverse lookup table
     for a directory of sift files."""
  cellid = getcellid(directory)
  data_out = getfile(directory, cellid + '-sift.npy')
  pydata_out = getfile(directory, cellid + '-pydata.npy')
  if not os.path.exists(data_out) or not os.path.exists(pydata_out):
    INFO('finding sift files')
    npy_save_sift_directory(directory, cellid)
  return np.load(data_out, mmap_mode=MMAP_MODE), np.load(pydata_out).item()

if __name__ == '__main__':
  from sys import argv
  data, mapping = npy_cached_load(argv[1])
  INFO('loaded %d features from %d images' % (len(data), len(mapping)))

# vim: et sw=2
