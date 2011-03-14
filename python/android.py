#!/usr/bin/env python
# reads xml and image files from the android Imageotag app
# (which includes tons of sensor data)

from config import *
import sys
import os
import xml.parsers.expat

class FlatXMLParser:
  """Just sets key=value for the innermost elements."""
  def __init__(self, xmlfile):
    self.current = None
    self.data = {}
    p = xml.parsers.expat.ParserCreate()
    p.StartElementHandler = self.start_element
    p.CharacterDataHandler = self.char_data
    p.ParseFile(open(xmlfile))

  def start_element(self, name, attrs):
    self.current = name

  def char_data(self, data):
    self.data[self.current] = data

class TaggedImage(object):
  def __init__(self, basedir, xmlfile, jpgfile):
    self.id = 'DSC_' + os.path.basename(xmlfile)[24:30].translate(None, '-_')
    self.basedir = basedir
    self.jpg = jpgfile
    self.sift = jpgfile[:-4] + 'sift.txt'
    xml = FlatXMLParser(xmlfile)
    self.data = xml.data
    self.time = int(self.data['clock_datetime'])
    self.gps_acc = float(self.data['gps_accuracy'])
    self.alt = self.gps_alt = float(self.data['gps_altitude'])
    self.lat = self.gps_lat = float(self.data['gps_latitude'])
    self.lon = self.gps_lon = float(self.data['gps_longitude'])
    self.net_acc = float(self.data['net_accuracy'])
    self.net_lat = float(self.data['net_latitude'])
    self.net_lon = float(self.data['net_longitude'])
    self.pitch = float(self.data['orientation_y'])
    self.yaw = float(self.data['orientation_x'])
    self.roll = float(self.data['orientation_z'])

  def getCanonicalName(self):
    """Name in line with previous naming schemes."""
    return "%s,%.6f,%.6f" % (self.id, self.gps_lat, self.gps_lon)

  def __str__(self):
    out = self.getCanonicalName()
    out += "\n file: %s" % self.jpg
    out += "\n time: %s" % self.time
    out += "\n gps accuracy: %s" % self.gps_acc
    out += "\n lat: %s" % self.gps_lat
    out += "\n lon: %s" % self.gps_lon
    out += "\n alt: %s" % self.gps_alt
    out += "\n net accuracy: %s" % self.net_acc
    out += "\n net lat: %s" % self.net_lat
    out += "\n net lon: %s" % self.net_lon
    out += "\n pitch: %s" % self.pitch
    out += "\n yaw: %s" % self.yaw
    out += "\n roll: %s" % self.roll
    return out

class AndroidReader:
  def __init__(self, basedir):
    """Reads image data from basedir"""
    self.basedir = basedir
    self.images = []
    xml_files = sorted(filter(lambda k: k.endswith('.xml'), os.listdir(basedir)))
    for xml in xml_files:
      name = xml[10:-4]
      jpg = name + '.jpg'
      if not os.path.exists(os.path.join(basedir, jpg)):
        INFO('W: could not find %s' % jpg)
        continue
      image = TaggedImage(basedir, os.path.join(basedir, xml), jpg)
      self.images.append(image)

  def __iter__(self):
    return self.images.__iter__()

  def __str__(self):
    mn = sum(map(lambda i: i.net_acc, self.images))/len(self.images)
    mg = sum(map(lambda i: i.gps_acc, self.images))/len(self.images)
    return "%d images, mean gps acc %f, net acc %f" % (len(self.images), mg, mn)

if __name__ == '__main__':
  if len(sys.argv) < 2 or not os.path.isdir(sys.argv[1]):
    print "Usage: %s <directory>" % sys.argv[0]
    exit()
  reader = AndroidReader(sys.argv[1])
  for image in reader:
    print image
  print reader

# vim: et sw=2
