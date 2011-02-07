import info
import math
import numpy as np
import numpy.linalg as linalg

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
    self.kv = tuple(sorted(kv.iteritems()))
    self.filteredlen = None

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
  def __init__(self, taglist):
    self.tags = []
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
