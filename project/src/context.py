# the Query class usually named Q
# holds all information relevant to a query image

# the Context class usually named C
# holds configuration information needed to run the query

import os
import info
import Image
import pixels
import numpy as np
import util
import shutil
import query
from config import *

from android import AndroidReader
from tags import TagCollection

class _Query:
  def __init__(self, instance=None):
    self.jpgpath = None
    self.siftpath = None
    self.sensor_lat = None
    self.sensor_lon = None
    self._query_lat = None
    self._query_lon = None
    self.pgm_scale = 1.0
    self.data = {}
    if instance:
      self.__dict__.update(instance.__dict__)

  def setQueryCoord(self, qlat, qlon):
    self._query_lat, self._query_lon = qlat, qlon

  def setSensorCoord(self, lat, lon):
    self.sensor_lat, self.sensor_lon = lat, lon

  def copy(self):
    return _Query(self)
  
  @property
  def query_lat(self):
    return self._query_lat or self.sensor_lat

  @property
  def query_lon(self):
    return self._query_lon or self.sensor_lon

  @property
  def name(self):
    return self.jpgname[:-4]
  
  @property
  def siftname(self):
    return os.path.basename(self.siftpath)

  @property
  def jpgname(self):
    return os.path.basename(self.jpgpath)
  
  def check(self):
    assert os.path.exists(self.jpgpath)
    assert os.path.exists(self.siftpath)
    assert self.sensor_lat is not None and self.sensor_lon is not None
    assert self.jpgname[:-4] == self.siftname[:-8]
    
class _Context(object):

  def __init__(self, context=None):
    # set frozen attr
    self.unfreeze()

    # default parameters
    self.QUERY = None # set before calling characterize()
    self.params = query.PARAMS_DEFAULT.copy()
    self.cacheEnable = 0 # instance-local caching of results
    self.do_posit = 0
    self.solve_pnp = 0
    self.print_per = 1
    self.amb_cutoff = None
    self.amb_padding = 50
    self.one_big_cell = 0
    self.added_error = 0
    self.restrict_cells = False
    self.override_cells = False
    self.max_matches_to_analyze = 1
    self.disable_filter_step = False
    self.stop_on_homTrue = 0
    self.put_into_dirs = 0
    self.locator_function = lambda C, Q: [(Q.sensor_lat, Q.sensor_lon)]
    self.weight_by_distance = False
    self.weight_by_coverage = False
    self._dbdir = False
    self.cellradius = 236.6
    self.shuffle_cells = False
    self.ranking_min_consistent = 10
    self.overlap_method = None
    self.ranking_max_considered = 100
    self.spatial_comb = 0
    self.match_callback = None
    self.dump_hom = 0
    self.solve_pose = 0
    self.log_failures = True
    self.solve_bad = 0
    self.ambiguity = 75
    self._test_r = None
    self._test_d = None
    self.datasource = None
    self.matchdistance = 25
    self.selection = None
    self.tagcompute = True # false is like NO_HOM, NO_DRAW
    self.show_feature_pairs = False
    self.compute2dpose = False
    self.min_reproj = False
    self.ncells = 10 # if ambiguity<100, 9 is max possible by geometry
    self.verbosity = 1
    self.resultsdir = os.path.expanduser('~/topmatches')
    self.reproj_file = os.path.expanduser('~/topmatches/reproj')
    self.topnresults = []
    self.maindir = os.path.expanduser('/media/DATAPART2')
    self.aarondir='fuzz2'
    self.bundler = False
    self.bundlerdir = '/media/DATAPART2/ah/bundler'

    # lazy load
    self._tags = None
    self._pixelmap = None
    
    # pull in new data
    if context:
      self.__dict__.update(context.__dict__)

    # always starts writable
    self.unfreeze()

  def __setattr__(self, *args):
    if self.frozen:
      raise Exception("Assignment to read-only context")
    object.__setattr__(self, *args)

  def freeze(self):
    object.__setattr__(self, 'frozen', True)

  def unfreeze(self):
    object.__setattr__(self, 'frozen', False)

  def frozen_copy(self):
    copy = self.copy()
    copy.freeze()
    return copy

  def copy(self):
    return _Context(self)

  def loadkey(self, key):
    p = os.path.join(CACHE_PATH, key + '.npy')
    if os.path.exists(p):
      return np.load(p).item()
    else:
      return None

  def savekey(self, key, value):
    p = os.path.join(CACHE_PATH, key + '.npy')
    save_atomic(lambda d: np.save(d, value), p)

  def pickleable(self):
    self.tags # force lazy load here
    copy = self.copy()
    for k,v in self.__dict__.items():
      if type(v) == type(lambda:0):
        copy.__setattr__(k, None)
    return copy

  @property
  def tags(self):
    if not self._tags:
      src = os.path.dirname(__file__)
      self._tags = TagCollection(os.path.join(src, 'tags-canonical.csv'))
#      src = os.path.dirname(os.path.dirname(__file__))
#      self._tags = TagCollection(
#        os.path.join(src, 'tags.csv'),
#        os.path.join(src, 'bearings.csv'),
#      )
    return self._tags

  @property
  def pixelmap(self):
    if not self._pixelmap:
      self._pixelmap = pixels.PixelMap(os.path.join(self.maindir, self.dbdump))
    return self._pixelmap

  @property
  def dbdump(self):
    if self.QUERY == 'emeryville':
      return os.path.join(self.maindir, 'Research/cells/emeryville/link_to_single_cell')
    elif self.QUERY == 'oakland1' or self.QUERY == 'oak-test':
      return '/media/DATAPART1/oakland/earthmine/rect' #os.path.join(self.maindir, 'Research/collected_images/earthmine-oakland/oakland-rect')
    elif self.QUERY == 'cory-4':
      return os.path.join(self.maindir, 'Research/collected_images/cory/db-4')
    elif self.QUERY == 'cory-25':
      return os.path.join(self.maindir, 'Research/collected_images/cory/db-25')
    elif self.QUERY == 'cory-2':
      return os.path.join(self.maindir, 'Research/collected_images/cory/db-2')
    elif self.QUERY == 'cory-5':
      return os.path.join(self.maindir, 'Research/collected_images/cory/db-5')
    return os.path.join(self.maindir, 'Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829')

  @property
  def dbdir(self):
    if self._dbdir:
      return self._dbdir
    elif self.QUERY == 'emeryville':
      return os.path.join(self.maindir, 'Research/cells/emeryville/single/')
    elif self.QUERY == 'oakland1':
      return '/media/DATAPART1/oakland/cells'
    elif self.QUERY == 'oak-test':
      if self._test_r == self._test_d == 236.6:
        return '/media/DATAPART1/oakland/cells'
      else:
        return '/media/DATAPART1/oak,r=%s,d=%s' % (self._test_r, self._test_d)
    elif self.QUERY == 'q5-test' or self.QUERY == 'q4-test':
      if self._test_r == self._test_d == 236.6:
        return '/media/DATAPART1/earthmine-fa10.1-culled,r=d=236.6'
      else:
        return '/media/DATAPART1/earthmine-fa10.1-culled,r=%s,d=%s' % (self._test_r, self._test_d)
    elif self.QUERY == 'cory-4':
      return    os.path.join(self.maindir, 'Research/cells/cory-4')
    elif self.QUERY == 'cory-25':
      return    os.path.join(self.maindir, 'Research/cells/cory-25')
    elif self.QUERY == 'cory-2':
      return    os.path.join(self.maindir, 'Research/cells/cory-2')
    elif self.QUERY == 'cory-5':
      return    os.path.join(self.maindir, 'Research/cells/cory-5')
    return os.path.join(self.maindir, 'Research/cells/g=100,r=d=236.6/')

  @property
  def fuzzydir(self):
    return os.path.join(self.maindir, 'fuzzylocs/%s' % self.QUERY)

  @property
  def matchdir(self):
    if self.QUERY == 'q5-test' or self.QUERY == 'q4-test' or self.QUERY == 'oak-test':
      celldesc = [y for y in self.dbdir.split('/') if y][-1]
      return os.path.join(self.maindir, 'Research/results/%s/matchescells(%s),%s,%s' % (self.QUERY, celldesc, self.QUERY, query.searchtype(self.params)))
    return os.path.join(self.maindir, 'Research/results/%s/matchescells(g=100,r=d=236.6),%s,%s' % (self.QUERY, self.QUERY, query.searchtype(self.params)))

  @property
  def infodir(self):
    if self.QUERY == 'oakland1' or self.QUERY == 'oak-test':
      return self.dbdump
    return os.path.join(self.maindir, 'Research/collected_images/earthmine-fa10.1/37.871955,-122.270829')

  @property
  def querydir(self):
    if self.QUERY == 'oakland1' or self.QUERY == 'oak-test':
      return '/media/DATAPART1/oakland/query/set1'
    return os.path.join(self.maindir, '%s/' % self.QUERY)

  @property
  def drawtopcorr(self):
    return self.tagcompute and 'NO_DRAW' not in os.environ

  @property
  def compute_hom(self):
    return self.tagcompute and 'NO_HOM' not in os.environ

  def initdirs(self):
    """Creates and cleans result data directories."""
    if not os.path.exists(self.matchdir):
        os.makedirs(self.matchdir)
    if self.drawtopcorr:
        if os.path.exists(self.resultsdir):
            shutil.rmtree(self.resultsdir)
        os.makedirs(self.resultsdir)
    if self.solve_pose:
        if os.path.exists(self.pose_param['resultsdir']):
            shutil.rmtree(self.pose_param['resultsdir'])
        os.makedirs(self.pose_param['resultsdir'])

  def iter_queries(self):
    for Q in self.iter_queries_unfiltered():
      if not self.selection:
        yield Q
      else:
        for s in self.selection:
          if s in Q.name:
            yield Q
            break

  def iter_queries_unfiltered(self):
    """Returns iter over _Query for files in query"""

    #if query taken from a cell phone
    if self.QUERY == 'query4' or self.QUERY == 'query4-cropped' or \
        self.QUERY == 'query4a' or self.QUERY == 'query5horizontal' or \
        self.QUERY == 'q5-test' or self.QUERY == 'q4-test' or \
        self.QUERY == 'query5vertical' or self.QUERY == 'oakland1' \
        or self.QUERY == 'oak-test':
      def iter0():
        for a in AndroidReader(self.querydir):
          image = _Query()
          image.siftpath = os.path.join(a.basedir, a.sift)
          image.jpgpath = os.path.join(a.basedir, a.jpg)
          image.setSensorCoord(
            *info.add_error((a.lat, a.lon), self.added_error)
          )
          if self.QUERY == 'query5horizontal' or self.QUERY == 'oakland1' or self.QUERY == 'oak-test':
            image.pgm_scale = 512/1952.0
          image.check()
          image.datasource = a
          yield image
      return iter0()

    #if taken from d2x with no coordinates
    if self.QUERY == 'emeryville' or self.QUERY == 'cory-25' or self.QUERY == 'cory-2' or self.QUERY == 'cory-5':
      def iter1():
        for file in util.getSiftFileNames(self.querydir):
          image = _Query()
          image.siftpath = os.path.join(self.querydir, file)
          #TODO: this uses a hard coded dimension
          image.pgm_scale = max(Image.open(os.path.join(self.querydir, image.siftname[:-8] + '.pgm')).size) / max(2592.0, 1456.0)
          image.jpgpath = os.path.join(self.querydir, image.siftname[:-8] + '.jpg')
          image.setSensorCoord(0,0)
          image.check()
          yield image
      return iter1()

    def iter2():
      for file in util.getSiftFileNames(self.querydir):
        image = _Query()
        if self.QUERY == 'query2':
          image.pgm_scale = 512/2592.0
        image.siftpath = os.path.join(self.querydir, file)
        image.jpgpath = os.path.join(self.querydir, image.siftname[:-8] + '.JPG')
        image.setSensorCoord(*info.add_error(info.getQuerySIFTCoord(file), self.added_error))
        if self.QUERY != 'query4-matlab':
          image.check()
        yield image
    return iter2()

DEFAULT_CONTEXT = _Context().frozen_copy()

# vim: et sw=2
