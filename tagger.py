from info import distance
import Image, ImageDraw
import json
import numpy.linalg as linalg
import numpy as np
import math

# all angles are as in unit circle; unit degrees

class TagCollection:
  def __init__(self, taglist):
    self.tags = []
    for line in open(taglist):
      line = line.split(',')
      tag = {
        'lat': line[0],
        'lon': line[1],
        'alt': line[2],
      }
      key = None
      for elem in line[3:]
        if key is None:
          key = elem
        else:
          tag[key] = elem
          key = None
    self.tags.append(tag)

  def select_frustrum(self, lat, lon, yaw, max_angle=90, max_radius=50):
    contained = []
    for tag in tags:
      lat2 = tag['lat']
      lon2 = tag['lon']
      dist = distance(lat, lon, lat2, lon2)
      if dist > max_radius:
        continue
      v1 = (math.cos(math.radians(yaw)), math.sin(math.radians(yaw)))
      v2 = (lon,lat)/linalg.norm((lon2-lon, lat2-lat))
      degrees = math.acos(np.dot(v1,v2))*180/math.pi
      if degrees > max_angle:
        continue
      contained.append(tag)
    return contained

class TaggedImage:
  def __init__(self, image, info, db):
    self.db = db
    with open(info) as f:
      self.info = json.JSONDecoder().decode(f.read())
    self.image = Image.open(image)

  def get_frustrum(self):
    return self.db.select_frustrum(self.info['lat'], self.info['lon'], self.info['yaw'])

  def get_distance(self, tag, pixel):
    conn = ddObject()
    viewId = self.info['id']
    viewPixels = [ddViewPixel(pixel[0], pixel[1]))]
    response = conn.getLocationsForViewPixels(viewID, viewPixels)
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
    tag_positions = {}
    for pixel in self.get_pixels():
      for tag in possible_tags:
        dist = self.get_distance(tag, pixel)
        if dist < DIST_THRESHOLD:
          tag_positions[tag] = (dist, pixel)
    for tag in tag_positions:
      places = tag_positions[tag]
      if places:
        yield (tag, min(places)[1])

  def save(self, output):
    points = self.get_tag_points()
    output = ImageDraw.Draw(self.image)
    for tag, coord in points:
      draw.text(coord, str(tag))
    self.image.save(output, 'jpg')

if __name__ == '__main__':
  pass # TODO test image

# vim: et sw=2
