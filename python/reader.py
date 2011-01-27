#!/usr/bin/env python
# Loads feature files, caching feature vectors for future runs.
# Use get_reader(type).load_file(file)
# Use get_reader(type).load_cell(directory)
#
#  typical numpy datatype stored: {
#    'names': ['vec', 'geom', 'index'],
#    'formats': ['128uint8', '4float32', 'uint16'],
#  }

from config import *
import numpy as np
import time
import os

def get_reader(typehint):
  if 'sift' in typehint:
    return SIFTReader()
  elif 'surf' in typehint:
    return SURFReader()
  elif 'chog' in typehint:
    return CHoGReader()
  else:
    raise NotImplementedError

class FeatureReader(object):
  def __init__(self, name, dtype):
    self.name = name
    self.dtype = dtype

  def is_feature_file(self, filename):
    raise NotImplementedError

  def write_features_to_array(self, filename, offset, dataset, index):
    raise NotImplementedError

  def count_features_in_file(self, filename):
    raise NotImplementedError

  def get_corresponding_jpg(self, filename):
    return filename.split(self.name)[0] + '.jpg'

  def get_feature_files_in_dir(self, directory):
    for name in os.listdir(directory):
      if self.is_feature_file(name):
        yield os.path.join(directory, name)

  def count_features_in_dir(self, directory):
    num_features = 0
    for f in self.get_feature_files_in_dir(directory):
      num_features += self.count_features_in_file(f)
    return num_features

  def save_directory(self, directory, cellid):
    """Writes all features found in a directory to a file.
       Also builds a reverse lookup table."""

    num_features = self.count_features_in_dir(directory)
    assert num_features > 0
    dataset = np.ndarray(num_features, self.dtype)
    offset = 0
    image_index = 0
    lookup_table = {}
    # now begin the actual read
    start = time.time()
    for name in os.listdir(directory):
      if self.is_feature_file(name):
        offset = self.write_features_to_array(os.path.join(directory, name),
          offset, dataset, image_index)
        lookup_table[image_index] = name
        if image_index % 200 == 0:
          INFO('%d/%d features read' % (offset, num_features))
          INFO_TIMING('%d features per second' % (offset/(time.time() - start)))
        image_index += 1
    INFO('saving features...')
    for dest in getdests(directory, cellid + ('-%s.npy' % self.name)):
      save_atomic(lambda d: np.save(d, dataset), dest)
    for dest in getdests(directory, cellid + ('-%s-pydata.npy' % self.name)):
      save_atomic(lambda d: np.save(d, lookup_table), dest)

  def load_file(self, file):
    """Reads features from file."""
    num_features = self.count_features_in_file(file)
    assert num_features > 0
    dataset = np.ndarray(num_features, self.dtype)
    self.write_features_to_array(file, 0, dataset, 0)
    return dataset

  def load_cell(self, directory):
    """Efficiently loads a matrix of features and reverse lookup table
       for a directory of files."""
    cellid = getcellid(directory)
    #features, geom, index
    data_out = getfile(directory, cellid + ('-%s.npy' % self.name))
    #mapping
    pydata_out = getfile(directory, cellid + ('-%s-pydata.npy' % self.name))
    if not os.path.exists(data_out) or not os.path.exists(pydata_out):
      INFO('finding %s files' % self.name)
      self.save_directory(directory, cellid)
    return np.load(data_out), np.load(pydata_out).item()

class SIFTReader(FeatureReader):
  sift_dtype = {
    'names': ['vec', 'geom', 'index'],
    'formats': ['128uint8', '4float32', 'uint16'],
  }

  def __init__(self):
    super(SIFTReader, self).__init__('sift', self.sift_dtype)

  def is_feature_file(self, filename):
    return 'sift.txt' in filename

  def count_features_in_file(self, siftfile):
    with open(siftfile) as f:
      return int(f.readline().split()[0])

  def sift_iterator(self, siftname):
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

  def write_features_to_array(self, siftname, offset, dataset, index):
    """Adds features from a sift file to a dataset.
       Returns the new offset into the matrix."""
    step = 0
    for chunk in self.sift_iterator(siftname):
      if chunk.dtype == np.float32:
        dataset['geom'][offset] = chunk
      else: # vector info
        dataset['vec'][offset][step:step+len(chunk)] = chunk
        step += len(chunk)
        if step >= 128:
          dataset['index'][offset] = index
          step = 0 # on to next vector
          offset += 1
    return offset

class CHoGReader(FeatureReader):
  chog_dtype = {
    'names': ['vec', 'geom', 'index'],
    'formats': ['63float', '5float32', 'uint16'],
  }

  def __init__(self):
    super(CHoGReader, self).__init__('chog', self.chog_dtype)

  def is_feature_file(self, filename):
    return 'chog.txt' in filename

  def count_features_in_file(self, chogfile):
    with open(chogfile) as f:
      return f.read().count('\n')

  def chog_iterator(self, chogname):
    with open(chogname) as data:
      for line in data:
        yield np.fromstring(line, sep=' ', dtype=np.float)

  def write_features_to_array(self, chogname, offset, dataset, index):
    for chunk in self.chog_iterator(chogname):
      assert len(chunk) == 68
      dataset['geom'][offset] = chunk[:5]
      dataset['vec'][offset] = chunk[5:]
      dataset['index'][offset] = index
      offset += 1
    return offset

class SURFReader(FeatureReader):
  surf_dtype = {
    'names': ['vec', 'geom', 'index'],
    'formats': ['64int32', '3uint16', 'uint16'],
  }

  def __init__(self):
    super(SURFReader, self).__init__('surf', self.surf_dtype)

  def is_feature_file(self, filename):
    return 'surf.npy' in filename

  def count_features_in_file(self, surffile):
    return len(np.load(surffile))

  def write_features_to_array(self, surfname, offset, dataset, index):
    data = np.load(surfname)
    # quantize...
    dataset[offset:offset+len(data)]['vec'] = np.int32(127*data['vec'])
    dataset[offset:offset+len(data)]['geom'] = data['geom']
    dataset[offset:offset+len(data)]['index'] = [index]*len(data)
    return offset + len(data)

# vim: et sw=2