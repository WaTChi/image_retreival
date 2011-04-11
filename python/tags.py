import info
import os
import math
import geom
import numpy as np
import numpy.linalg as linalg
from earthMine import ddObject, ddGetViews

class OcclusionSummary:
  def __init__(self, ocs):
    self.ocs = ocs

  def bestView(self, viewId, lat, lon, distance_f):
    for o in self.ocs:
      if o['id'] == viewId:
        print "--> EXACT MATCH"
        return o
    p = []
    for o in self.ocs:
      dist = distance_f(o['loc'])
      p.append((dist, o))
    p.sort()
    print "--> LAT LON DISTANCE MATCH", p[0][0]
    return p[0][1]
  
  def hasNonOccludedView(self, viewId, lat, lon, distance_f, thresh):
    p = []
    for o in self.ocs:
      dist = distance_f(o['loc'])
      p.append((dist, o))
    p.sort()
    for d, o in p:
      if d > thresh:
        break
      if not o['occ']:
        return True
    return False

class Tag:
  """Representation of an EarthMine tag."""
  def __init__(self, kv):
    self.lat = kv['lat']
    self.lon = kv['lon']
    self.alt = kv['alt']
    del kv['lat']
    del kv['lon']
    del kv['alt']
    if kv['name'] == 'business':
      self.business = True
      del kv['name']
    else:
      self.business = False
    self.name = kv.get('Name')
    if not self.name:
      self.name = kv.get('name')
    self.kv = tuple(sorted(kv.iteritems()))
    self.filteredlen = None
    self._ocs = None

  def getOcclusionSummary(self, C):
    if self._ocs:
      return self._ocs
    radius = 100
    maxResults = 100
    key = "%f,%f,%d,%d.views" % (self.lat, self.lon, radius, maxResults)
    val = C.loadkey(key)
    if val:
      self._ocs = val
    else:
      conn = ddObject()
      result = ddGetViews(conn, self.lat, self.lon, radius=radius, maxResults=maxResults)
      summary = []
      for i, r in enumerate(result):
        item = {
                  'loc': r['view-location'],
                  'dir': r['view-direction'],
                  'occ': r['is-known-occluded'],
                  'id': r['id']
               }
        summary.append(item)
      self._ocs = OcclusionSummary(summary)
      C.savekey(key, self._ocs)
    return self._ocs

  # info = EarthMineImageInfo
  def emIsVisible(self, source, C, thresh):
    ocs = self.getOcclusionSummary(C)
    def distance_f(oloc):
      return info.distance(source.lat, source.lon, oloc['lat'], oloc['lon'])
    return ocs.hasNonOccludedView(source.viewId, source.lat, source.lon, distance_f, thresh)

  def isVisible2(self, source, tree2d, elat, elon):
    numoccs = self.howOccluded(source, tree2d, elat, elon)
    return numoccs <= 0

  # info = EarthmineImageInfo
  # tree3d defines method:
  #     countHigherPtsNear(lat, lon, alt, threshold) -> int
  def howOccluded(self, source, tree3d, elat, elon):

    # subdivide line until granularity of 1 meter is reached
    def recur(lat1, lon1, alt1, lat2, lon2, alt2):
      mlat, mlon, malt = geom.midpoint(lat1, lon1, alt1, lat2, lon2, alt2)
      if geom.distance3d6(self.lat, self.lon, self.alt, mlat, mlon, malt) < 4.0:
        return 0
      elif geom.distance3d6(lat1, lon1, alt1, lat2, lon2, alt2) < 1.0:
        return tree3d.countHigherPtsNear(mlat, mlon, malt, 2.0)
      else:
        return recur(lat1, lon1, alt1, mlat, mlon, malt) +\
               recur(lat2, lon2, alt2, mlat, mlon, malt)

    return recur(self.lat, self.lon, self.alt, elat, elon, source.alt)

  def xydistance(self, d2):
    x1, y1, z1 = self.lat, self.lon, self.alt
    x2, y2, z2 = d2['lat'], d2['lon'], d2['alt']
    xydist = info.distance(x1, y1, x2, y2)
    return xydist

  def distance(self, d2):
    xydist = self.xydistance(d2)
    vert = abs(self.alt-d2['alt'])
    return math.sqrt(xydist**2 + vert**2)

  def __str__(self):
    return str(self.kv)

  def __repr__(self):
    return str(self)

  def __hash__(self):
    return hash(self.lat) + hash(self.lon) + hash(self.alt) + hash(self.business) + hash(self.kv)

  def __equals__(self, other):
    return self.lat == other.lat and self.lon == other.lon and self.alt == other.alt and self.business == other.business and self.kv == other.kv

  def filteritems(self):
    out = {}
    kv = dict(self.kv)
    for k,v in kv.iteritems():
      if kv[k]:
        out[k] = kv[k]
    return out

  def __iter__(self):
    ordered = ['name', 'Name', 'House Number', 'Phone Number', 'Comment', 'Field 1', 'Field 2', 'Field 3', 'Field 4']
    tag = self.filteritems()
    for k in ordered:
      if k in tag:
        if type(k) == str:
          yield str(tag[k])
    for k in tag:
      if k not in ordered:
        yield str(tag[k])

  def __len__(self):
    if self.filteredlen is None:
      self.filteredlen = len(self.filteritems())
    return self.filteredlen

class TagCollection:
  """Parses EarthMine's tag export format."""

  def googlefmt(self):
    with open("out.googleearth", 'w') as file:
      for tag in self.tags:
        print >>file, ",%s,%f,%f" % (tag.name, tag.lat, tag.lon)

  def __init__(self, taglist):
    self.tags = []
    if not os.path.exists(taglist):
      print "W: tag file not found!"
      return
    for line in open(taglist):
      line = line.split(',')
      tag = {
        'lat': float(line[0]),
        'lon': float(line[1]),
        'alt': float(line[2]),
        'name': line[3],
      }
      key = None
      for elem in line[4:]:
        if key is None:
          key = elem
        else:
          tag[key] = elem.strip()
          key = None
      self.tags.append(Tag(tag))

  def select_frustum(self, lat, lon, yaw, fov=270, radius=100):
    contained = []
    for tag in self.tags:
      dist = info.distance(lat, lon, tag.lat, tag.lon)
      if dist > radius:
        continue
      v1 = (math.sin(math.radians(yaw)), math.cos(math.radians(yaw)))
      v2 = (tag.lon-lon, tag.lat-lat)/linalg.norm((tag.lon-lon, tag.lat-lat))
      degrees = math.acos(np.dot(v1,v2))*180/math.pi
      assert degrees > 0
      if degrees > fov/2:
        continue
      contained.append(tag)
    return contained

# vim: et sw=2
