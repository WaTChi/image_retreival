import os
import info
import Image
import pixels
import util
import shutil
import query

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
    self.ransac_min_filt = 1
    self.do_posit = 0
    self.print_per = 1
    self.max_matches_to_analyze = 1
    self.stop_on_homTrue = 0
    self.put_into_dirs = 0
    self.locator_function = lambda C, Q: [(Q.sensor_lat, Q.sensor_lon)]
    self.cellradius = 236.6
    self.match_callback = None
    self.dump_hom = 0
    self.ambiguity = 75
    self.datasource = None
    self.matchdistance = 25
    self.ncells = 10 # if ambiguity<100, 9 is max possible by geometry
    self.verbosity = 1
    self.resultsdir = os.path.expanduser('~/topmatches')
    self.topnresults = []
    self.maindir = os.path.expanduser('/media/DATAPART2')

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

  def pickleable(self):
    copy = self.copy()
    for k,v in self.__dict__.items():
      if type(v) == type(lambda:0):
        copy.__setattr__(k, None)
    return copy

  @property
  def tags(self):
    if not self._tags:
      self._tags = TagCollection(os.path.join(self.maindir, 'Research/app/code/tags.csv'))
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
    return os.path.join(self.maindir, 'Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829')

  @property
  def dbdir(self):
    if self.QUERY == 'emeryville':
      return os.path.join(self.maindir, 'Research/cells/emeryville/single/')
    return os.path.join(self.maindir, 'Research/cells/g=100,r=d=236.6/')

  @property
  def fuzzydir(self):
    return os.path.join(self.maindir, 'fuzzylocs/%s' % self.QUERY)

  @property
  def matchdir(self):
    return os.path.join(self.maindir, 'Research/results/%s/matchescells(g=100,r=d=236.6),%s,%s' % (self.QUERY, self.QUERY, query.searchtype(self.params)))

  @property
  def infodir(self):
	return os.path.join(self.maindir, 'Research/collected_images/earthmine-fa10.1/37.871955,-122.270829')

  @property
  def querydir(self):
    return os.path.join(self.maindir, '%s/' % self.QUERY)

  @property
  def drawtopcorr(self):
    return 'NO_DRAW' not in os.environ

  @property
  def compute_hom(self):
    return 'NO_HOM' not in os.environ

  def initdirs(self):
    """Creates and cleans result data directories."""
    if not os.path.exists(self.matchdir):
        os.makedirs(self.matchdir)
    if self.drawtopcorr:
        if os.path.exists(self.resultsdir):
            shutil.rmtree(self.resultsdir)
        os.makedirs(self.resultsdir)

  def iter_queries(self):
    """Returns iter over _Query for files in query"""
    if self.QUERY == 'query4':
      def iter0():
        for a in AndroidReader(self.querydir):
          image = _Query()
          image.siftpath = os.path.join(a.basedir, a.sift)
          image.jpgpath = os.path.join(a.basedir, a.jpg)
          image.setSensorCoord(a.lat, a.lon)
          image.check()
          image.datasource = a
          yield image
      return iter0()

    if self.QUERY == 'emeryville':
      def iter1():
        for file in util.getSiftFileNames(self.querydir):
          image = _Query()
          image.siftpath = os.path.join(self.querydir, file)
          image.pgm_scale = max(Image.open(os.path.join(self.querydir, image.siftname[:-8] + '.pgm')).size) / max(2592.0, 1456.0)
          image.jpgpath = os.path.join(self.querydir, image.siftname[:-8] + '.jpg')
          image.setSensorCoord(0,0)
          image.check()
          yield image
      return iter1()

    def iter2():
      for file in util.getSiftFileNames(self.querydir):
        image = _Query()
        image.siftpath = os.path.join(self.querydir, file)
        image.jpgpath = os.path.join(self.querydir, image.siftname[:-8] + '.JPG')
        image.setSensorCoord(*info.getQuerySIFTCoord(file))
        if self.QUERY != 'query4-matlab':
          image.check()
        yield image
    return iter2()

DEFAULT_CONTEXT = _Context().frozen_copy()

# vim: et sw=2
