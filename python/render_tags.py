#!/usr/bin/env python
# all angles are as in unit circle; unit degrees

from config import INFO
import Image, ImageDraw, ImageFont
import os
import earthMine
import geom
import reader
import util
import colorsys
import math
import info
import numpy as np
from tags import TagCollection
from android import AndroidReader

def get_image_info(db_img):
  if os.path.basename(db_img).startswith('DSC_'):
    return NikonImageInfo(db_img)
  return EarthmineImageInfo(db_img, db_img[:-4] + '.info')

class ImageInfo:
  """Metadata source for image. Provides:
      self.image self.focal_length
      self.lat self.lon self.alt
      self.pitch self.yaw self.roll"""

  def get_loc_dict(self):
    return {'lat': self.lat, 'lon': self.lon, 'alt': self.alt}

  def get_pixel_location(self, pixel):
    out = self.get_pixel_locations([pixel])
    if out is None:
      return None
    return out.get(pixel, None)

  def get_pixel_locations(self, pixels):
    """Return None if not supported."""
    raise NotImplementedError

class ComputedImageInfo(ImageInfo):
  """computes pose from db matches"""
  def __init__(self, image, lat, lon):
    self.image = Image.open(image)
    self.lat = lat
    self.lon = lon

  def get_pixel_locations(self, pixels):
    return None

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

class NikonImageInfo(ImageInfo):
  def __init__(self, image):
    self.lat = 0
    self.lon = 0
    self.alt = 0
    self.yaw = 0
    self.roll = 0
    self.viewId = 0
    self.image = Image.open(image)
    self.focal_length = 700 # XXX arbitrary!
    center = geom.center((0,0), self.image.size)
    self.pitch = 0

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
    self.roll = 0
    self.viewId = self.info['id']
    self.image = Image.open(image) if image else None
    if not self.image:
      return
    self.focal_length = 0.8625 * self.image.size[0]
    center = geom.center((0,0), self.image.size)
    cloc = self.get_pixel_location(center)
    if not cloc:
      self.pitch = 0
      INFO("WARNING: Earthmine returned None for center pixel; set pitch=0")
      return
    hyp = geom.distance3d(cloc, self.info['view-location'])
    d_alt = cloc['alt'] - self.alt
    self.pitch = np.arcsin(d_alt/hyp)

  def get_pixel_locations(self, pixels):
    return earthMine.ddGetAllPixels(pixels, self.viewId)
         
class TaggedImage:
  """Tags an EarthMine image."""
  def __init__(self, image, source, db):
    self.db = db
    self.source = source
    self.focal_length = 700 # overwritten by source update below
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

  def map_tags_camera(self, elat=0, elon=0, ep=0, ey=0, er=0):
    "Returns (tag, (dist, pixel)) pairs using camera transform."
    if not elat or not elon:
      elat = self.lat
      elon = self.lon
    tags = []
    possible_tags = self.get_frustum()

    for tag in possible_tags:
      pz, px = geom.lltom(elat, elon, tag.lat, tag.lon)
      py = tag.alt - self.alt;
      x, y, z = geom.camera_transform(px, py, pz, self.pitch + ep, self.yaw + ey, self.roll + er)
      coord = geom.project2d(x, y, z, self.source.focal_length)
      pixel = geom.center(coord, self.image.size)
#      tags.append((tag, (0, geom.constrain(pixel, self.image.size))))
      tags.append((tag, (0, pixel)))

    return tags

  def map_tags_ocs(self, C):

    tags = self.map_tags_camera()
    accepted, bad = [], []

    for (tag, (_, pixel)) in tags:
      if tag.emIsVisible(self.source, C, 20):
        accepted.append((tag, (_, pixel)))
      else:
        bad.append((tag, (999, pixel)))

    return accepted + bad

  def map_tags_hybrid(self, pixelmap, C, elat, elon):
    """Uses tag projection from estimated lat, lon.
       Tags are filtered by using the image's pt cloud
       when source tag is visible in the db image. Otherwise,
       3d occlusion detection is performed."""
    THRESHOLD = 15.0 # meters
    tags = self.map_tags_camera(elat, elon)
    accepted = []
    outside = []
    bad = []
    for (tag, (_, pixel)) in tags:
      location = pixelmap[geom.picknearest(pixelmap, *pixel)]
      if location is None:
        if not geom.contains(pixel, self.image.size):
          outside.append((tag, (_, pixel)))
        else:
          bad.append((tag, (999, pixel)))
      else:
        dist = tag.xydistance(location)
        if dist < THRESHOLD:
          accepted.append((tag, (_, pixel)))
        elif not geom.contains(pixel, self.image.size):
          outside.append((tag, (_, pixel)))
        else:
          bad.append((tag, (999, pixel)))

    # use ocs method for tags outside db image
    cell = util.getclosestcell(self.lat, self.lon, C.dbdir)[0]
    cellpath = os.path.join(C.dbdir, cell)
    tree3d = reader.get_reader(C.params['descriptor']).load_tree3d(cellpath, C)
    for (tag, (_, pixel)) in outside:
      if tag.isVisible2(self.source, tree3d, elat, elon):
        accepted.append((tag, (_, pixel)))
      else:
        bad.append((tag, (15, pixel)))

    return accepted + bad

  def map_tags_hybrid2(self, pixelmap, C):
    """Tags are filtered by using the image's pt cloud
       when source tag is visible in the db image. Otherwise,
       a combination of earthmine occlusion queries and
       database occlusion queries are performed."""
    tags = self.map_tags_camera(self.lat, self.lon)
    accepted = []
    outside = []
    bad = []
    THRESHOLD = 15
    outside = tags # XXX
#    for (tag, (_, pixel)) in tags:
#      location = pixelmap[geom.picknearest(pixelmap, *pixel)]
#      if location is None:
#        if not geom.contains(pixel, self.image.size):
#          outside.append((tag, (_, pixel)))
#        else:
#          bad.append((tag, (999, pixel)))
#      else:
#        dist = tag.xydistance(location)
#        if dist < THRESHOLD:
#          accepted.append((tag, (_, pixel)))
#        elif not geom.contains(pixel, self.image.size):
#          outside.append((tag, (_, pixel)))
#        else:
#          bad.append((tag, (999, pixel)))

    cell = util.getclosestcell(self.lat, self.lon, C.dbdir)[0]
    cellpath = os.path.join(C.dbdir, cell)
    pm = reader.get_reader(C.params['descriptor'])\
      .load_PointToViewsMap(cellpath, C.infodir)

    for (tag, (_, pixel)) in outside:
      vis, t = pm.hasView(C, tag.lat, tag.lon, self.lat, self.lon, self.yaw, 30)
      emv = tag.emIsVisible(self.source, C, 30)
      if (vis or emv):
        if geom.norm_compatible(tag, self):
          accepted.append((tag, (_, pixel)))
        else:
          bad.append((tag, (12, pixel)))
      else:
        bad.append((tag, (17, pixel)))

    return accepted + bad

  def map_tags_lookup(self, C):
    tags = self.map_tags_camera()
    cell = util.getclosestcell(self.lat, self.lon, C.dbdir)[0]
    cellpath = os.path.join(C.dbdir, cell)

    pm = reader.get_reader(C.params['descriptor'])\
      .load_PointToViewsMap(cellpath, C.infodir)

    accepted, bad = [], []

    for (tag, (_, pixel)) in tags:
      vis, t = pm.hasView(C, tag.lat, tag.lon, self.lat, self.lon, self.yaw, 20)
      if vis:
        accepted.append((tag, (_, pixel)))
      else:
        bad.append((tag, (999, pixel)))
    return accepted + bad

  def map_tags_hybrid3(self, pixelmap, C):
    THRESHOLD = 15.0 # meters
    tags = self.map_tags_camera()
    accepted = []
    outside = []
    bad = []
    min_upper_bound = 0
    obs = self.source.get_loc_dict()
    for (tag, (_, pixel)) in tags:
      location = pixelmap[geom.picknearest(pixelmap, *pixel)]
      if location is None:
        if not geom.contains(pixel, self.image.size):
          outside.append((tag, (_, pixel)))
        else:
          bad.append((tag, (999, pixel)))
      else:
        dist = tag.xydistance(location)
        if dist < THRESHOLD:
          accepted.append((tag, (_, pixel)))
          min_upper_bound = max(min_upper_bound, tag.distance(obs))
        elif not geom.contains(pixel, self.image.size):
          outside.append((tag, (_, pixel)))
        else:
          bad.append((tag, (999, pixel)))

    cell = util.getclosestcell(self.lat, self.lon, C.dbdir)[0]
    cellpath = os.path.join(C.dbdir, cell)
    pm = reader.get_reader(C.params['descriptor'])\
      .load_PointToViewsMap(cellpath, C.infodir)

    for (tag, (_, pixel)) in outside:
      vis, t = pm.hasView(C, tag.lat, tag.lon,\
        self.lat, self.lon, self.yaw, 30)
      emv = tag.emIsVisible(self.source, C, 30)
      bv = tag.distance(obs) < min_upper_bound + THRESHOLD
      if (vis or emv or bv):
        if geom.norm_compatible(tag, self):
          accepted.append((tag, (_, pixel)))
        else:
          bad.append((tag, (12, pixel)))
      else:
        accepted.append((tag, (15, pixel)))
    return accepted + bad

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

  def taggedcopy(self, points, image):
    MIN_SIZE = 1
    draw = ImageDraw.Draw(image)
    def dist(p):
      return info.distance(p[0].lat, p[0].lon, self.lat, self.lon)
    points.sort(key=dist, reverse=True) # draw distant first
    for tag, (dist, point) in points:
      color = self.colordist(dist, 30.0)
      size = int(300.0/info.distance(tag.lat, tag.lon, self.lat, self.lon))
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
      img = image.copy()
      draw2 = ImageDraw.Draw(img)
      draw2.rectangle([top_left, bottom_right], fill='#000')

      image = Image.blend(image, img, 0.75)
      draw = ImageDraw.Draw(image)
      # end black container
      draw.ellipse((point[0]-size/2,point[1]-size/2,point[0]+size/2,point[1]+size/2), fill=color)
      for line in tag:
        draw.text((point[0] + off_x, point[1] + off_y), line, fill=color, font=font)
        off_y += max(size, MIN_SIZE)
#      if dist:
#        INFO('mapping tag at %f meters error' % dist)
    return image

  def draw(self, points, output):
    image = self.taggedcopy(points, self.image)
    image.save(output, 'png')
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
  jpg = 'testdata/x8.jpg'
  source = EarthmineImageInfo(jpg, jpg[:-4] + '.info')
  img = TaggedImage(jpg, source, db)
  points = img.map_tags_camera()
  img.draw(points, 'testdata/out.png')

if __name__ == '__main__':
  _test2()

# vim: et sw=2
