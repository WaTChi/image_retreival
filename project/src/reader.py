#!/usr/bin/env python
# Loads feature files, caching feature vectors for future runs.
# Use get_reader(type).load_file(file)
# Use get_reader(type).load_cell(directory)
#
#  typical numpy datatype stored: {
#    'names': ['vec', 'geom', 'index'],
#    'formats': ['128uint8', '4float32', 'uint16'],
#  }

import random
import itertools
from config import *
import query
import numpy as np
import shutil
import pyflann
import raster_tools
import time
from info import distance
import geom
import render_tags
import sys
import os

# memory hungry! beware.
memoize = False
ans = {}

def ichunks(iterable, size):
    it = iter(iterable)
    i = 0
    chunk = list(itertools.islice(it,size))
    while chunk:
        yield (i, chunk)
        i += size
        chunk = list(itertools.islice(it,size))

def config_mem_pin(val):
  global memoize
  memoize = val

def hashcells(cellpath):
  if type(cellpath) is list:
    return tuple(sorted(cellpath))
  return cellpath

def get_reader(typehint):
  if 'sift' in typehint:
    return SIFTReader()
  elif 'surf' in typehint:
    return SURFReader()
  elif 'chog' in typehint:
    return CHoGReader()
  else:
    raise NotImplementedError

LOCATION_UNIT = 1/1e5 # 1 meter
disc = lambda la: int(la/LOCATION_UNIT)

# provides map from 2dpt -> earthmine views of point
class PointToViewsMap(object):

  def __init__(self, lookup_table):
    self.lookup_table = lookup_table
    INFO_TIMING("initialized pt to view map with %d entries" % len(lookup_table))

  def _expandSet(self, k, di):
    for dx in range(-di, di + 1):
      for dy in range(-di, di + 1):
        yield (k[0] + dx, k[1] + dy)

  def _getViewStrs(self, lat, lon, rad=5):
    k0 = disc(lat), disc(lon)
    imgset = set()
    for k in self._expandSet(k0, rad):
      imgset = imgset.union(self.lookup_table.get(k, set()))
    return imgset

  def getViews(self, C, lat, lon, rad):
    for dbimg in self._getViewStrs(lat, lon, rad):
      info = os.path.join(C.infodir, os.path.basename(dbimg)[:-8] + '.info')
      source = render_tags.EarthmineImageInfo(None, info)
      yield source

  # return True if (lat, lon) is visible from (qlat, qlon, qyaw)
  # yaw is in radians
  def hasView(self, C, lat, lon, qlat, qlon, qyaw, thresh):
    dist = distance(lat, lon, qlat, qlon)
    rad = min(20, max(3, int(dist/10)))

    def similar(view):
      return distance(view.lat, view.lon, qlat, qlon) < thresh

    def yawof(view):
      return view.yaw

    views = list(self.getViews(C, lat, lon, rad))
#    print "dist is %d, rad is %d, nviews is %d" % (dist, rad, len(views))
    closeyaws = map(lambda (k,v): v, sorted([(distance(v.lat, v.lon, qlat, qlon), yawof(v)) for v in views]))[:10]
    norm = geom.compute_norm(closeyaws)

    return any([similar(v) for v in views]), norm

class FeatureReader(object):
  def __init__(self, name, dtype):
    self.name = name
    self.dtype = dtype
    self.filters = {
      'locally_significant': self.filter_locally_significant,
      'locally_insignificant': self.filter_locally_insignificant,
      'random_half': self.filter_random_half,
      'random_quarter': self.filter_random_quarter,
    }

  def filter_locally_insignificant(self, fgen, d0, m0, q):
    return self.filter_locally_significant(fgen, d0, m0, q, True)

  def filter_locally_significant(self, fgen, d0, m0, q, rev=False):
    T = 70000 # pts ABOVE threshold allowed
    R = 25.0 # max visibility distance
    NN = 10 # neighbors to consider

    flann = pyflann.FLANN()
    INFO("building flann index")
    q.flann_setup_index(flann, d0, m0, None)
    INFO("building 3d map")
    map3d = self.load_3dmap_for_cell(q.cellpath, d0, m0, q.infodir)
    INFO("begin filtering features")
    params = q.params.copy()
    params['num_neighbors'] = NN
    for offset, fchunk in fgen:
      fchunk_filtered = []
      arr = np.ndarray(len(fchunk), self.dtype)
      arr[:] = fchunk
      results, dists = flann.nn_index(arr['vec'], **params)
      # for each feature
      for _, (index_arr, dist_arr) in enumerate(zip(results, dists)):
        refpt = map3d[offset + _]
        # for each match of feature
        has_r = False
        ok = False
        for i, d in zip(index_arr, dist_arr):
          pt = map3d[i]
          r = distance(pt['lat'], pt['lon'], refpt['lat'], refpt['lon'])
          if r > R:
            has_r = True
            if (not rev and d > T) or (rev and d < T):
              ok = True
              break
        if ok or not has_r:
          fchunk_filtered.append(fchunk[_])
      print 'filter chunk %d/%d' % (len(fchunk_filtered), len(fchunk))
      yield fchunk_filtered

  def filter_random_quarter(self, fgen, d0, m0, q):
    for offset, fchunk in fgen:
      arr = []
      for i, x in enumerate(fchunk):
        if i % 4 == 0:
          arr.append(x)
      yield arr

  def filter_random_half(self, fgen, d0, m0, q):
    for offset, fchunk in fgen:
      arr = []
      for i, x in enumerate(fchunk):
        if i % 2 == 0:
          arr.append(x)
      yield arr

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

  def count_features_in_files(self, files):
    num_features = 0
    for f in files:
      num_features += self.count_features_in_file(f)
    return num_features

  def save_shuffled(self, dbdir, dirs, supercell, overlap_method):
    """Shuffles files in cells and places them in supercell"""
    files = []
    closed = set()
    for dir in [os.path.join(dbdir, d) for d in dirs]:
      for file in os.listdir(dir):
        if self.is_feature_file(file):
          if overlap_method is None or file not in closed:
            closed.add(file)
            files.append(os.path.join(dir, file))
    random.shuffle(files)

    # list<set<name>, set<path>> indexed by cell
    sets = []
    for i in range(len(dirs)):
      sets.append((set(), set()))

    # divide them up into cells
    while files:
      f = files.pop()
      name = os.path.basename(f)
      while True:
        names, paths = random.choice(sets)
        if name not in names:
          names.add(name)
          paths.add(f)
          break
      # supersample an average of `overlap` times
      if overlap_method:
        extra_samples = max(0, random.gauss(overlap_method, 0.2))
        for i in range(int(round(extra_samples))):
          iter = 0
          while True:
            if iter > 50:
              INFO("W: ITER EXCEEDED BOUNDS")
              break
            iter += 1
            names, paths = random.choice(sets)
            if name not in names:
              names.add(name)
              paths.add(f)
              break

    for i, old_name in enumerate(map(os.path.basename, dirs)):
      new_cell = os.path.join(supercell, old_name)
      if os.path.exists(new_cell):
        shutil.rmtree(new_cell)
      os.mkdir(new_cell)
      print "Building new cell of size %d" % len(sets[i][0])
      for file in sets[i][1]:
        os.symlink(file,
          os.path.join(new_cell, os.path.basename(file)))

  def save_multi(self, dirs, cellid):
    """Generalization of the old save_directory.
       Writes all features found in dirs to a file.
       Deduplicate by files.
       Also builds a reverse lookup table."""

    files = set()
    closed = set()
    import config
    for dir in dirs:
      for file in os.listdir(dir):
        if self.is_feature_file(file):
          if config.using_weighted_union() or file not in closed:
            files.add(os.path.join(dir, file))
            closed.add(file)

    num_features = self.count_features_in_files(files)
    assert num_features > 0
    dataset = np.ndarray(num_features, self.dtype)
    offset = 0
    image_index = 0
    lookup_table = {}
    # now begin the actual read
    start = time.time()
    for file in files:
      name = os.path.basename(file)
      offset = self.write_features_to_array(file, offset,
        dataset, image_index)
      lookup_table[image_index] = name
      if image_index % 200 == 0:
        INFO('%d/%d features read' % (offset, num_features))
        INFO_TIMING('%d features per second' % (offset/(time.time() - start)))
      image_index += 1
    INFO('saving features... [%s]' % cellid)
    for dest in getdests(dirs[0], cellid + ('-%s.npy' % self.name)):
      save_atomic(lambda d: np.save(d, dataset), dest)
    for dest in getdests(dirs[0], cellid + ('-%s-pydata.npy' % self.name)):
      save_atomic(lambda d: np.save(d, lookup_table), dest)

  def load_tree3d(self, directory, C):
    """Returns a tree class against which queries of the form
       tree3d.countPtsNear(lat, lon, threshold_meters)
       can be performed efficiently."""
    pixmap_dir = C.infodir
    try:
      map3d = self.load_3dmap_for_cell(directory, None, None, None)
    except:
      INFO('exception... rereading 3d pts')
      dataset, mapping = self.load_cell(directory)
      INFO('building 3d map...')
      map3d = self.load_3dmap_for_cell(directory, dataset, mapping, pixmap_dir)

    amap3d = map3d.view((np.float64, 3))
    return query.Tree3D(amap3d, C, getcellid(directory))

  # these maps are on the order of 5MB in size
  # and contain ~100k entries
  def load_PointToViewsMap(self, directory, pxdir):

    cellid = getcellid(directory)
    v_out = getfile(directory, cellid + '-pvmap.npy')

    if os.path.exists(v_out):
      return PointToViewsMap(np.load(v_out).item())

    dataset, mapping = self.load_cell(directory)
    map3d = self.load_3dmap_for_cell(directory, dataset, mapping, pxdir)

    lookup_table = {}
    for i in range(len(map3d)):
      if i % 123456 == 0:
        print ('\r -> %d/%d' % (i, len(map3d))),
        sys.stdout.flush()
      ds_lat = disc(map3d[i]['lat'])
      ds_lon = disc(map3d[i]['lon'])
      key = (ds_lat, ds_lon)
      if key not in lookup_table:
        lookup_table[key] = set()
      lookup_table[key].add(mapping[dataset[i]['index']])

    for dest in getdests(directory, cellid + '-pvmap.npy'):
      save_atomic(lambda d: np.save(d, lookup_table), dest)

    return PointToViewsMap(lookup_table)

  def load_hsv_for_cell(self, directory, dataset, mapping, pixmap_dir):
    """For the cell specified, return a vector v such that
       hsv,depth,plane(dataset[i]) == v[i]
       This vector (~30MB) will be cached on disk."""
    hsv_dtype = {
      'names': ['h3', 's3', 'v3',
                'plane', 'depth'],
      'formats': ['uint8', 'uint8', 'uint8',
                  'uint32', 'float32'],
    }
    cellid = getcellid(directory)
    v_out = getfile(directory, cellid + '-hsv.npy')
    if os.path.exists(v_out):
      return np.load(v_out)
    print v_out, "does not exist"
    INFO('reading hsv, depth, plane values for %s' % cellid)

    hsv = np.zeros(len(dataset), hsv_dtype)
    oldindex = -1
    for i, row in enumerate(dataset):
      coord = row['geom']
      index = row['index']
      if index != oldindex:
        if index % 200 == 0:
          INFO('read %d/%d pts' % (i, len(dataset)))
        oldindex = index
        base = mapping[index][:-len('sift.txt')]
        raster = os.path.join(pixmap_dir, base + '.jpg')
        plane = os.path.join(pixmap_dir, base + '-planes.png')
        depth = os.path.join(pixmap_dir, base + '-depth.png')
        raster = raster_tools.RasterImage(raster)
        if os.path.exists(plane):
          plane = raster_tools.RasterImage(plane)
        else:
          plane = None
        if os.path.exists(depth):
          depth = raster_tools.RasterImage(depth)
        else:
          depth = None
      for X in [3]:
        h,s,v = raster.sampleHSV(coord[0], coord[1], X)
        hsv[i]['h%d' % X] = h
        hsv[i]['s%d' % X] = s
        hsv[i]['v%d' % X] = v
      if plane:
        hsv[i]['depth'] = depth.getDepth(coord[0], coord[1])
      if depth:
        hsv[i]['plane'] = plane.getPlane(coord[0], coord[1])
    for dest in getdests(directory, cellid + '-hsv.npy'):
      save_atomic(lambda d: np.save(d, hsv), dest)
    return hsv

  def load_3dmap_for_cell(self, directory, dataset, mapping, pixmap_dir):
    """For the cell specified, return a vector v such that
       3dcoord(dataset[i]) == v[i]
       This vector (~100MB) will be cached on disk."""
    # IMPORTANT:
    # IF ANY FIELD IS SET TO ZERO, THE DISTANCE IS UNKNOWN/FAR AWAY
    map3d_dtype = {
      'names': ['lat', 'lon', 'alt'],
      'formats': ['float64', 'float64', 'float64'],
    }
    cellid = getcellid(directory)
    v_out = getfile(directory, cellid + '-map3d.npy')
    if os.path.exists(v_out):
      return np.load(v_out)
    print v_out, "does not exist"
    INFO('reading 3d pts for %s' % cellid)

    import pixels # avoid circular imports

    pixmap = pixels.PixelMap(pixmap_dir)
    map3d = np.zeros(len(dataset), map3d_dtype)
    oldindex = -1
    for i, row in enumerate(dataset):
      coord = row['geom']
      index = row['index']
      if index != oldindex:
        if index % 200 == 0:
          INFO('read %d/%d pts' % (i, len(dataset)))
        oldindex = index
        pts3d = pixmap.open(mapping[index])
      pt3d = pts3d[int(coord[0]), int(coord[1])]
      if pt3d is not None:
        map3d[i]['lat'] = pt3d['lat']
        map3d[i]['lon'] = pt3d['lon']
        map3d[i]['alt'] = pt3d['alt']
    for dest in getdests(directory, cellid + '-map3d.npy'):
      save_atomic(lambda d: np.save(d, map3d), dest)
    return map3d

  def load_file(self, file):
    """Reads features from file."""
    num_features = self.count_features_in_file(file)
    assert num_features > 0
    dataset = np.ndarray(num_features, self.dtype)
    self.write_features_to_array(file, 0, dataset, 0)
    return dataset

  def get_supercelldir(self, dbdir, directories, overlap_method):
    """return directory of shuffled cells"""
    if overlap_method is not None:
      suffix = '_overlap_%f' % overlap_method
    else:
      suffix = ''
    supercell = '/media/DATAPART1/supercells/' + getcellid(directories) + suffix
    stamp = os.path.join(supercell, 'stamp')
    if not os.path.exists(stamp):
      if not os.path.exists(supercell):
        os.mkdir(supercell)
      self.save_shuffled(dbdir, directories, supercell, overlap_method)
      open(stamp, 'w')
    return supercell

  def load_cell(self, directory, q, criteria=None):
    if criteria is None:
      return self._load_cell(directory)
    elif criteria in self.filters:
      d, m = self._load_cell(directory)
      d2, m2 = self.make_compressed_cell(d, m, criteria, directory, q)
      return d2, m2
    else:
      assert False, criteria

  def make_compressed_cell(self, d0, m0, criteria, directory, q):
    cellid = getcellid(directory)
    data_out = getfile(directory, cellid + ('-%s-%s-test2.npy' % (self.name, criteria)))
    if not os.path.exists(data_out):
      INFO("building compressed cell of locally significant features")
      tmp = np.ndarray(len(d0), self.dtype)
      offset = 0
      for fchunk in self.filters[criteria](ichunks(d0, 10000), d0, m0, q):
        tmp[offset:offset+len(fchunk)] = fchunk
        offset += len(fchunk)
      print "retained %d/%d features" % (offset, len(d0))
      # truncate tmp
      tmp = tmp[:offset]
      print "writing to", data_out
      save_atomic(lambda d: np.save(d, tmp), data_out)
      print "write done"

    # index mapping does not change
    d2 = np.load(data_out)
    assert len(d2) > 0, data_out
    return d2, m0

  def _load_cell(self, directory):
    """Efficiently loads a matrix of features and reverse lookup table
       for a directory of files."""
    global ans
    key = hashcells(directory)
    if key in ans:
      INFO('*** MEMOIZED CELL ***')
      return ans[key]
    cellid = getcellid(directory)
    #features, geom, index
    data_out = getfile(directory, cellid + ('-%s.npy' % self.name))
    #mapping
    pydata_out = getfile(directory, cellid + ('-%s-pydata.npy' % self.name))
    if not os.path.exists(data_out) or not os.path.exists(pydata_out):
      INFO('finding %s files' % self.name)
      if type(directory) is list:
        self.save_multi(directory, cellid)
      else:
        self.save_multi([directory], cellid)
    d, m = np.load(data_out), np.load(pydata_out).item()
    if memoize:
      ans[key] = (d,m)
    return d, m

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
    print siftfile
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
