#!/usr/bin/env python
# all angles are as in unit circle; unit degrees

from config import INFO
import Image, ImageDraw, ImageFont
import earthMine
import geom
import colorsys
import math
import info
import numpy as np
from tags import TagCollection
from android import AndroidReader

class ImageInfo:
  """Metadata source for image. Provides:
      self.image self.focal_length
      self.lat self.lon self.alt
      self.pitch self.yaw self.roll"""

  def get_pixel_location(self, pixel):
    out = self.get_pixel_locations([pixel])
    if out is None:
      return None
    return out.get(pixel, None)

  def get_pixel_locations(self, pixels):
    """Return None if not supported."""
    assert False, "NotImplemented"

class AndroidImageInfo(ImageInfo):
  """for images from AndroidReader"""
  def __init__(self, image):
    self.__dict__.update(image.__dict__)
    self.image = Image.open(image.jpg)
    # TODO not sure if these are right.
    # It's hard to tell when lat/lon/alt are so wrong
    self.pitch = (self.pitch + 90) * math.pi / 180
    self.yaw = (self.yaw - 60) * math.pi / 180
    self.roll = (self.roll) * math.pi / 180
  
  def get_pixel_locations(self, pixels):
    return None

class EarthmineImageInfo(ImageInfo):
  def __init__(self, image, infofile):
    with open(infofile) as f:
      self.info = eval(f.read())
    self.lat = self.info['view-location']['lat']
    self.lon = self.info['view-location']['lon']
    self.alt = self.info['view-location']['alt']
    self.yaw = self.info['view-direction']['yaw']
    self.yaw = -(self.yaw+180)*math.pi/180
    self.roll = 0
    self.viewId = self.info['id']
    self.image = Image.open(image)
    self.focal_length = 0.8625 * self.image.size[0]
    center = geom.center((0,0), self.image.size)
    cloc = self.get_pixel_location(center)
    hyp = geom.distance3d(cloc, self.info['view-location'])
    d_alt = cloc['alt'] - self.alt
    self.pitch = np.arcsin(d_alt/hyp)
  
  def get_pixel_locations(self, pixels):
    """fetches an arbitrary amount of pixels from EarthMine Direct Data"""
    conn = earthMine.ddObject()
    viewPixels = [earthMine.ddViewPixel(p[0], p[1]) for p in pixels]
    locs = {}
    while viewPixels:
      response = conn.getLocationsForViewPixels(self.viewId, viewPixels[:490])
      viewPixels = viewPixels[490:] # limit for api
      for pixel, loc in response:
        if loc: # has valid mapping
          locs[(pixel['x'], pixel['y'])] = loc
    return locs # map of (x,y) to coords

class TaggedImage:
  """Tags an EarthMine image."""
  def __init__(self, image, source, db):
    self.db = db
    self.source = source
    self.focal_length = 700
    self.__dict__.update(self.source.__dict__)

  def get_frustum(self):
    return self.db.select_frustum(self.lat, self.lon, self.yaw)

  def save_point_cloud(self, output):
    cloud = self.source.get_pixel_locations(1,1)
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

  def map_tags_camera(self):
    "Returns (tag, (dist, pixel)) pairs using camera transform."

    THRESHOLD = 20.0 # meters
    tags = []
    possible_tags = self.get_frustum()

    for tag in possible_tags:
      pz, px = geom.lltom(self.lat, self.lon, tag.lat, tag.lon)
      py = tag.alt - self.alt;
      x, y, z = geom.camera_transform(px, py, pz, self.pitch, self.yaw, self.roll)
      coord = geom.project2d(x, y, z, self.focal_length)
      pixel = geom.center(coord, self.image.size)
      tags.append((tag, (0, geom.constrain(pixel, self.image.size))))

    # add some distance debugging info from earthmine
    debugtags = []
    locs = self.source.get_pixel_locations(map(lambda t: t[1][1], tags))
    if locs is None: # does not support
      return tags

    for (tag, (nulldist, pixel)) in tags:
      dist = tag.xydistance(locs.get(pixel, {'lat': 9999, 'lon': 9999, 'alt': 9999}))
      if dist < THRESHOLD:
        debugtags.append((tag, (dist, pixel)))
    return debugtags

  def map_tags_earthmine(self):
    "Returns (tag, (dist, pixel)) pairs using earthmine pixel data."
    THRESHOLD = 10.0
    possible_tags = self.get_frustum()
    locs = self.source.get_pixel_locations(list(self.get_pixels()))
    assert locs is not None
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
        tags.append((tag, mapping[tag]))
    INFO("mapped %d/%d possible tags" % (len(tags), len(possible_tags)))
    return tags

  def draw(self, points, output):
    draw = ImageDraw.Draw(self.image)
    MIN_SIZE = 1
    points.sort(key=lambda p: p[1][0], reverse=True) # draw distant first
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
      if dist:
        INFO('mapping tag at %f meters error' % dist)
    self.image.save(output, 'png')
    INFO("saved to %s" % output)

def _test():
  db = TagCollection('testdata/tags.csv')
  import os
  idir = 'testdata/input'
  odir = 'testdata/output'
  distlog = []
  for f in os.listdir(idir):
    if '.jpg' in f:
      output = os.path.join(odir, f[:-4] + '.png')
      jpg = os.path.join(idir, f)
      source = EarthmineImageInfo(jpg, jpg[:-4] + '.info')
      img = TaggedImage(jpg, source, db)
      points = img.map_tags_camera()
      for tag, (dist, point) in points:
        distlog.append(dist)
      img.draw(points, output)
  for i in [1, 10, 50, 100, len(distlog)]:
    INFO('top %d error is %f at %d samples' % (i, sum(sorted(distlog)[:i])/i, len(distlog)))

def _testq4():
  db = TagCollection('testdata/tags.csv')
  import os
  idir = 'testdata/query4'
  odir = 'testdata/q4_output'
  reader = AndroidReader(idir)
  for image in reader.images:
    output = os.path.join(odir, image.getCanonicalName().replace('.jpg', '.png'))
    source = AndroidImageInfo(image)
    img = TaggedImage(image.jpg, source, db)
    points = img.map_tags_camera()
    img.draw(points, output)

def _test2():
  db = TagCollection('testdata/tags.csv')
  jpg = 'testdata/input/37.8719888495,-122.269869743-0002.jpg'
  source = EarthmineImageInfo(jpg, jpg[:-4] + '.info')
  img = TaggedImage(jpg, source, db)
  points = img.map_tags_camera()
  img.draw(points, 'testdata/out.png')

if __name__ == '__main__':
  _testq4()

# vim: et sw=2
