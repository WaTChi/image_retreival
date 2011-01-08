#!/usr/bin/env python
# all angles are as in unit circle; unit degrees

from info import distance
from config import *
import earthMine
import Image, ImageDraw
import numpy.linalg as linalg
import numpy as np
import math

class TagCollection:
  def __init__(self, taglist):
    self.tags = []

    def frozen(dictionary):
      return tuple(sorted(dictionary.iteritems()))

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
      self.tags.append(frozen(tag))

  def select_frustrum(self, lat, lon, yaw, fov=60, radius=50):
    contained = []
    for tag in self.tags:
      t = dict(tag)
      lat2 = t['lat']
      lon2 = t['lon']
      dist = distance(lat, lon, lat2, lon2)
      if dist > radius:
        continue
      v1 = (math.sin(math.radians(yaw)), math.cos(math.radians(yaw)))
      v2 = (lon2-lon, lat2-lat)/linalg.norm((lon2-lon, lat2-lat))
      degrees = math.acos(np.dot(v1,v2))*180/math.pi
      assert degrees > 0
      if degrees > fov/2:
        continue
      contained.append(tag)
    return contained

class TaggedImage:
  def __init__(self, image, info, db):
    self.db = db
    with open(info) as f:
      self.info = eval(f.read())
    self.image = Image.open(image)

  def get_frustrum(self):
    return self.db.select_frustrum(self.info['view-location']['lat'], self.info['view-location']['lon'], self.info['view-direction']['yaw'])

  def get_distance(self, tag, pixel):
    return 0 # XXX

  def get_distances(self, pixels):
    assert len(pixels) < 500 # earthmine max
    conn = earthMine.ddObject()
    viewId = self.info['id']
    viewPixels = [earthMine.ddViewPixel(p[0], p[1]) for p in pixels]
    response = conn.getLocationsForViewPixels(viewId, viewPixels)
    locs = {}
    for pixel, loc in response:
      if loc: # has valid mapping
        locs[(pixel['x'], pixel['y'])] = loc
    return locs

  def get_pixels(self):
    for x in range(0, self.image.size[0], 25):
      for y in range(0, self.image.size[1], 35):
        yield (x,y)

  def partial_manhattan_distance(self, d1, d2):
    x1, y1, z1 = d1['lat'], d1['lon'], d1['alt']
    x2, y2, z2 = d2['lat'], d2['lon'], d2['alt']
    xydist = distance(x1, y1, x2, y2)
    vert = abs(z1-z2)
    return math.sqrt(xydist**2 + vert**2)

  def get_tag_points(self):
    "Returns collection of (tag, pixel) pairs"
    possible_tags = self.get_frustrum()
    tag_positions = dict([(tag, (999999, None)) for tag in possible_tags])
    locs = self.get_distances(list(self.get_pixels()))
    for pixel in self.get_pixels():
      if pixel not in locs:
        continue
      for tag in possible_tags:
        dist = self.partial_manhattan_distance(locs[pixel], dict(tag))
        tag_positions[tag] = min(tag_positions[tag], (dist, pixel))
    tags = []
    for tag in tag_positions:
      places = tag_positions[tag]
      if places[1]:
        t = dict(tag)
        t['map_dist'], t['mapped_coord'] = places
        tags.append(t)
    INFO("mapped %d/%d possible tags" % (len(tags), len(possible_tags)))
    return tags

  def draw(self, points, output):
    draw = ImageDraw.Draw(self.image)
    for tag in points:
      if tag['map_dist'] < 10.0:
        coord = tag['mapped_coord']
        offset = -30
        for key in tag:
          draw.text((coord[0], coord[1] + offset), str(tag[key]))
          offset += 10
        INFO('mapping tag at %f meters error' % tag['map_dist'])
      else:
        INFO('throwing away tag %s' % str(tag))
    self.image.save(output, 'jpeg')
    INFO("saved to %s" % output)

def _test():
  name = 'x5'
  db = TagCollection('testdata/tags.csv')
  img = TaggedImage('testdata/%s.jpg' % name, 'testdata/%s.info' % name, db)
  points = img.get_tag_points()
  img.draw(points, 'testdata/output.jpg')

if __name__ == '__main__':
  _test()

# vim: et sw=2
