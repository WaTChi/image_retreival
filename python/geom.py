# place for geometric utility functions
from info import distance
import numpy as np
import math
from numpy import matrix, sin, cos, sqrt

def lltom(lat1, lon1, lat2, lon2):
  """Returns lat2,lon2 relative to lat1,lon1 in meters.
     Does not handle wraparound at 360."""
  if lat1 == lat2 and lon1 == lon2:
    return 0
  dx = distance(lat1, lon1, lat2, lon1) * (-1 if lat2 < lat1 else 1)
  dy = distance(lat1, lon1, lat1, lon2) * (-1 if lon2 < lon1 else 1)
  return [dx,dy]

def camera_transform(px, py, pz, pitch, yaw, roll):
  """Translates units from px, py, pz to camera perspective.
     Units are in radians and meters."""
  yaw = -(yaw+180)*math.pi/180
  p = np.matrix([px, py, pz]).transpose()
  x, y, z = pitch, yaw, roll
  A = matrix(((1, 0, 0), (0, cos(x), -sin(x)), (0, sin(x), cos(x))))
  B = matrix(((cos(y), 0, sin(y)), (0, 1, 0), (-sin(y), 0, cos(y))))
  C = matrix(((cos(z), -sin(z), 0), (sin(z), cos(z), 0), (0, 0, 1)))
  d = A*B*C*p
  return d[0].item(), d[1].item(), d[2].item()

def project2d(x, y, z, f):
  """Projects 3d coordinate onto 2d plane with focal length f."""
  x, y = x*f/z, y*f/z
  return x, y

def center(point, screen):
  """Translates point so that (0,0) is in center of screen."""
  return int(screen[0]/2 + point[0]), int(screen[1]/2 + point[1])

def constrain(point, screen):
  """Translates point so that it is within the screen."""
  return max(0, min(screen[0], point[0])), max(0, min(screen[1], point[1]))

def distance3d(d1, d2):
  """Finds distance between {'lat': deg, 'lon': deg, 'alt': meters} points."""
  x1, y1, z1 = d1['lat'], d1['lon'], d1['alt']
  x2, y2, z2 = d2['lat'], d2['lon'], d2['alt']
  xydist = distance(x1, y1, x2, y2)
  vert = abs(d1['alt']-d2['alt'])
  return sqrt(xydist**2 + vert**2)

# vim: et sw=2
