import os
import info
import util
import shutil
import query

from android import AndroidReader

class _Query:
  def __init__(self, instance=None):
    self.jpgpath = None
    self.siftpath = None
    self.sensor_lat = None
    self.sensor_lon = None
    self._query_lat = None
    self._query_lon = None
    self.data = {}
    if instance:
      self.__dict__.update(instance.__dict__)

  def setQueryCoord(self, qlat, qlon):
    self._query_lat, self._query_lon = qlat, qlon

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
    assert self.jpgname[:-4] == self.siftname[:-8]
    return self.jpgname[:-4]
  
  @property
  def siftname(self):
    return os.path.basename(self.siftpath)

  @property
  def jpgname(self):
    return os.path.basename(self.jpgpath)

class _Context(dict):

  def __init__(self, mapping):
    # default parameters
    self.QUERY = None # set before calling characterize()
    self.params = query.PARAMS_DEFAULT.copy()
    self.cacheEnable = 0 # instance-local caching of results
    self.ransac_min_filt = 1
    self.do_posit = 0
    self.print_per = 1
    self.num_images_to_print = 1
    self.corrfilter_printed = 0 # keep trying for homTrue until num_images_to_print
    self.put_into_dirs = 0
    self.showHom = 0
    self.locator_function = lambda image: [(image.lat, image.lon)]
    self.cellradius = 236.6
    self.match_callback = None
    self.ambiguity = 75
    self.matchdistance = 25
    self.ncells = 10 # if ambiguity<100, 9 is max possible by geometry
    self.verbosity = 1
    self.resultsdir = os.path.expanduser('~/topmatches')
    self.topnresults = []
    self.maindir = os.path.expanduser('/media/DATAPART2')

    # pull in new data
    self.__dict__.update(mapping)

    # always starts writable
    self.unfreeze()

  def __setattr__(self, *args):
    if self.frozen:
      raise Exception("Assignment to read-only context")
    dict.__setattr__(self, *args)

  def freeze(self):
    dict.__setattr__(self, 'frozen', True)

  def unfreeze(self):
    dict.__setattr__(self, 'frozen', False)

  def frozen_copy(self):
    copy = _Context(self)
    copy.freeze()
    return copy

  def copy(self):
    return _Context(self)

  @property
  def dbdump(self):
    if self.QUERY == 'emeryville':
      return os.path.join(self.maindir, 'Research/cells/g=100,r=d=236.6/0,0')
    return os.path.join(self.maindir, 'Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829')

  @property
  def fuzzydir(self):
    return os.path.join(self.maindir, 'fuzzylocs/%s' % self.QUERY)

  @property
  def matchdir(self):
    return os.path.join(self.maindir, 'Research/results/%s/matchescells(g=100,r=d=236.6),%s,%s' % (self.QUERY, self.QUERY, query.searchtype(self.params)))

  @property
  def querydir(self):
    return os.path.join(self.maindir, '%s/' % self.QUERY)

  def initdirs(self):
    """Creates and cleans result data directories."""
    if not os.path.exists(self.matchdir):
        os.makedirs(self.matchdir)
    if self.drawtopcorr:
        if os.path.exists(self.resultsdir):
            shutil.rmtree(self.resultsdir)
        os.makedirs(self.resultsdir)

  def iter_queries(self, querydir):
    """Returns iter over _Query for files in query"""
    if self.QUERY == 'query4':
      return AndroidReader(querydir)
    if self.QUERY == 'emeryville':
      def iter():
        for file in util.getSiftFileNames(self.querydir):
          image = _Query()
          image.siftpath = os.path.join(self.querydir, file)
          yield image
      return iter()
    def iter2():
      for file in util.getSiftFileNames(self.querydir):
        image = _Query()
        image.siftpath = os.path.join(self.querydir, file)
        image.setQueryCoord(info.getQuerySIFTCoord(file))
        yield image
    return iter2()

DEFAULT_CONTEXT = _Context().frozen_copy()

# vim: et sw=2
