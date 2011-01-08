#!/usr/bin/env python
# all angles are as in unit circle; unit degrees

from info import distance
import earthMine
import Image, ImageDraw
import numpy.linalg as linalg
import numpy as np
import math

def frozen(dictionary):
  return tuple(sorted(dictionary.iteritems()))

class TagCollection:
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
    conn = earthMine.ddObject()
    viewId = self.info['id']
    viewPixels = [earthMine.ddViewPixel(pixel[0], pixel[1])]
    response = conn.getLocationsForViewPixels(viewId, viewPixels)
    loc = response[pixel]
    return distance(loc[0], loc[1], tag['lat'], tag['lon'])

  def get_pixels(self):
    for x in range(self.image.size[0]):
      for y in range(self.image.size[1]):
        yield (x,y)

  def get_tag_points(self):
    "Returns collection of (tag, pixel) pairs"
    DIST_THRESHOLD = 1.0 # meters
    possible_tags = self.get_frustrum()
    tag_positions = dict([(tag, []) for tag in possible_tags])
    for pixel in self.get_pixels():
      for tag in possible_tags:
        dist = self.get_distance(tag, pixel)
        if dist < DIST_THRESHOLD:
          tag_positions[tag].append((dist, pixel))
    for tag in tag_positions:
      places = tag_positions[tag]
      if places:
        t = dict(tag)
        t['mapped_coord'] = min(places)[1]
        yield t

  def save(self, output):
    points = self.get_tag_points()
    draw = ImageDraw.Draw(self.image)
    for tag, coord in points:
      draw.text(coord, str(tag))
    self.image.save(output, 'jpg')

if __name__ == '__main__':
  db = TagCollection('tags.csv')
  img = TaggedImage('test.jpg', 'test.info', db)
  for tag in img.get_frustrum():
    print dict(tag)
  print list(img.get_tag_points())

# vim: et sw=2
