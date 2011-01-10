#!/usr/bin/env python
# all angles are as in unit circle; unit degrees

from config import *
import Image, ImageDraw, ImageFont
import colorsys
import earthMine
import geom
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
    xydist = self.xydistance(x1, y1, x2, y2)
    vert = abs(z1-z2)
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

  def select_frustum(self, lat, lon, yaw, fov=90, radius=100):
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

class TaggedImage:
  """Tags an EarthMine image."""
  def __init__(self, image, info, db):
    """db is a TagCollection
       image an EarthMine jpg
       info an EarthMine .info file"""
    self.db = db
    with open(info) as f:
      self.info = eval(f.read())
    self.lat = self.info['view-location']['lat']
    self.lon = self.info['view-location']['lon']
    self.alt = self.info['view-location']['alt']
    self.yaw = self.info['view-direction']['yaw']
    self.pitch = self.info['view-direction']['pitch']
    self.roll = 0
    self.field_of_view = self.info['field-of-view']
    self.viewId = self.info['id']
    self.image = Image.open(image)

  def get_frustum(self):
    return self.db.select_frustum(self.lat, self.lon, self.yaw)

  def get_pixel_locations(self, pixels):
    """fetches an arbitrary amount of pixels from EarthMine Direct Data"""
    INFO('fetching %d pixels' % len(pixels))
    conn = earthMine.ddObject()
    viewPixels = [earthMine.ddViewPixel(p[0], p[1]) for p in pixels]
    locs = {}
    while viewPixels:
      response = conn.getLocationsForViewPixels(self.viewId, viewPixels[:490])
      viewPixels = viewPixels[490:] # limit for api
      INFO('%d pixels more' % len(viewPixels))
      for pixel, loc in response:
        if loc: # has valid mapping
          locs[(pixel['x'], pixel['y'])] = loc
    return locs # map of (x,y) to coords

  def save_point_cloud(self, output):
    cloud = self.get_pixel_locations(1,1)
    np.save(output, cloud)

  def colordist(self, dist, bound=100.0):
    offset = .33 # green ... red
    return tuple(map(lambda i: int(i*255), colorsys.hsv_to_rgb(min(1, dist/bound + offset), 1, 1)))

  def colorcycle(self, dist, bound=10.0):
    return tuple(map(lambda i: int(i*255), colorsys.hsv_to_rgb(dist/bound % 1.0, 1, 1)))

  def visualize_point_cloud(self, cloudfile, output):
    cloud = np.load(cloudfile).item()
    img = self.image.copy()
    pixels = img.load()
    for x,y in cloud:
      d = cloud[x,y]
      lat, lon = d['lat'], d['lon']
      pixels[x,y] = self.colordist(info.distance(lat, lon, self.lat, self.lon))
    img = Image.blend(img, self.image, 0.5)
    img.save(output, 'png')

  def get_pixels(self, stepx=25, stepy=35):
    for x in range(0, self.image.size[0], stepx):
      for y in range(0, self.image.size[1], stepy):
        yield (x,y)

  def refine_point(self, tag, place):
    pixel = place[1]
    worst = place[0]
    buf = []
    for dx in range(-50,50,6):
      for dy in range(-50,50,6):
        buf.append((pixel[0] + dx, pixel[1] + dy))
    locs = self.get_pixel_locations(buf)
    for pixel in buf:
      if pixel not in locs:
        continue
      dist = tag.distance(locs[pixel])
      place = min(place, (dist, pixel))
    INFO("improved dist by %f (%d%%)" % (worst - place[0], 100*(worst-place[0])/worst))
    return place

  def map_tags_camera(self):
    "Returns (tag, (dist, pixel)) pairs using camera transform."

    THRESHOLD = 40.0 # meters
    tags = []
    possible_tags = self.get_frustum()

    for tag in possible_tags:
      pz, px = geom.lltom(self.lat, self.lon, tag.lat, tag.lon)
      py = tag.alt - self.alt - 2.0 # XXX adjust for height of sensor
      focal_length = 0.9 * self.image.size[0] # TODO measure
      x, y, z = geom.camera_transform(px, py, pz, self.pitch, -(self.yaw+180)*math.pi/180, self.roll)
      coord = geom.project2d(x, y, z, focal_length)
      pixel = geom.center(coord, self.image.size)
      tags.append((tag, (0, geom.constrain(pixel, self.image.size))))

    # add some distance debugging info from earthmine
    debugtags = []
    locs = self.get_pixel_locations(map(lambda t: t[1][1], tags))
    for (tag, (nulldist, pixel)) in tags:
      dist = tag.xydistance(locs.get(pixel, {'lat': 9999, 'lon': 9999, 'alt': 9999}))
      if dist < THRESHOLD:
        debugtags.append((tag, (dist, pixel)))
    return debugtags

  def map_tags_earthmine(self):
    "Returns (tag, (dist, pixel)) pairs using earthmine pixel data."
    THRESHOLD = 10.0
    possible_tags = self.get_frustum()
    locs = self.get_pixel_locations(list(self.get_pixels()))
    mapping = dict([(tag, (999999, None)) for tag in possible_tags])
    for pixel in self.get_pixels():
      if pixel not in locs:
        continue
      for tag in possible_tags:
        dist = tag.distance(locs[pixel])
        if dist < mapping[tag][0]:
          mapping[tag] = (dist, pixel)
    tags = []
    for tag in possible_tags:
      if mapping[tag][0] < THRESHOLD:
#        mapping[tag] = self.refine_point(tag)
        tags.append((tag, mapping[tag]))
    INFO("mapped %d/%d possible tags" % (len(tags), len(possible_tags)))
    return tags

  def draw(self, points, output):
    draw = ImageDraw.Draw(self.image)
    MIN_SIZE = 2
    for tag, (dist, point) in points:
      color = self.colordist(dist, 10.0)
      size = int(200.0/info.distance(tag.lat, tag.lon, self.lat, self.lon))
      fontPath = "/usr/share/fonts/truetype/ttf-dejavu/DejaVuSans-Bold.ttf"
      font = ImageFont.truetype(fontPath, max(size, MIN_SIZE))
      off_x = -size*2
      off_y = -size*(len(tag)+1)
      # start black container
      top_left = (point[0] + off_x - 3, point[1] + off_y - 1)
      w = 10
      for line in tag:
        w = max(w, draw.textsize(line, font)[0])
      bottom_right = (point[0] + off_x + w + 3, point[1] + off_y + max(size, MIN_SIZE)*len(tag) + 3)
      img = self.image.copy()
      draw2 = ImageDraw.Draw(img)
      draw2.rectangle([top_left, bottom_right], fill='#000')
      draw2.ellipse((point[0]-size,point[1]-size,point[0]+size,point[1]+size), fill='#000')
      self.image = Image.blend(self.image, img, 0.75)
      draw = ImageDraw.Draw(self.image)
      # end black container
      draw.ellipse((point[0]-size/2,point[1]-size/2,point[0]+size/2,point[1]+size/2), fill=color)
      for line in tag:
        draw.text((point[0] + off_x, point[1] + off_y), line, fill=color, font=font)
        off_y += max(size, MIN_SIZE)
      INFO('mapping tag at %f meters error' % dist)
    self.image.save(output, 'png')
    INFO("saved to %s" % output)

def _test():
  db = TagCollection('testdata/tags.csv')
  import os
  idir = 'testdata/input'
  odir = 'testdata/output'
  for f in os.listdir(idir):
    if '.jpg' in f:
      jpg = os.path.join(idir, f)
      img = TaggedImage(jpg, jpg[:-4] + '.info', db)
      points = img.map_tags_camera()
      img.draw(points, os.path.join(odir, f[:-4] + '.png'))

def _test2():
  db = TagCollection('testdata/tags.csv')
  name = 'a'
  jpg = 'testdata/%s.jpg' % name
  img = TaggedImage(jpg, jpg[:-4] + '.info', db)
  points = img.map_tags_camera()
  img.draw(points, 'testdata/%s-out.png' % name)

if __name__ == '__main__':
  _test()

# vim: et sw=2
