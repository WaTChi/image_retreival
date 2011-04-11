# place for geometric utility functions
from info import distance
import numpy as np
import scipy.optimize as sio
import math
from numpy import matrix, sin, cos, sqrt

def euclideandist(x,y,x2,y2):
  return ((x-x2)**2 + (y-y2)**2)**.5

# all units are in radians
# finds minimum angle difference between t1,t2
def anglediff(t1, t2):
  a = math.e ** (1.0j * t1)
  b = math.e ** (1.0j * t2)
  c = b/a
  return abs(np.arctan2(c.imag, c.real))

# yaws = angles looking torwards the point
# returns the angle a such that the sum of cubes of
# anglediff(a, yaw_i) is minimized
def compute_norm(yaws):
  def f(t):
    error = 0.0
    for y in yaws:
        error += anglediff(y,t)**3
    return error
  return sio.brute(f, [(0,2*math.pi)])[0]

# yaw: yaw of potential viewer
# return True if (vlat, vlon) is visible at an angle yaw
def norm_compatible(norm, yaw):
  return anglediff(norm, yaw)*180/math.pi < 90

def picknearest(dict2d, x, y):
  """Returns closest key in dict of (x,y) to (x,y)"""
  return min(map(lambda k: (euclideandist(k[0],k[1],x,y),k), dict2d))[1]

def midpoint(lat1, lon1, alt1, lat2, lon2, alt2):
  return (lat1 + lat2) / 2.0, (lon1 + lon2) / 2.0, (alt1 + alt2) / 2.0

def lltom(lat1, lon1, lat2, lon2):
  """Returns lat2,lon2 relative to lat1,lon1 in meters.
     Does not handle wraparound at 360."""
  if lat1 == lat2 and lon1 == lon2:
    return [0, 0]
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

def contains(point, screen):
  return point == constrain(point, screen)

def distance3d6(x1, y1, z1, x2, y2, z2):
  xydist = distance(x1, y1, x2, y2)
  vert = abs(z1-z2)
  return sqrt(xydist**2 + vert**2)

def distance3d(d1, d2):
  """Finds distance between {'lat': deg, 'lon': deg, 'alt': meters} points."""
  x1, y1, z1 = d1['lat'], d1['lon'], d1['alt']
  x2, y2, z2 = d2['lat'], d2['lon'], d2['alt']
  return distance3d6(x1, y1, z1, z2, y2, z2)

# vim: et sw=2
