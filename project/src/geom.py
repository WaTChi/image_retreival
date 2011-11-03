# place for geometric utility functions
# usually relating to lat/lon/angles/pixels
# also see util.py and info.py

from info import distance, getbearing
import numpy as np
import scipy.optimize as sio
import math
from numpy import matrix, sin, cos, sqrt, pi

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

# return True if pt is visible at from viewpt
def norm_compatible(pt, viewpt, verbose=False):

  # no bearing: corner, middle of street, etc
  if pt.bearing is None:
    return True

  # really close up!
  if distance(pt.lat, pt.lon, viewpt.lat, viewpt.lon) < 20:
    return True

  yaw = getbearing(pt.lat, pt.lon, viewpt.lat, viewpt.lon)
  diff = anglediff(pt.bearing*pi/180, yaw*pi/180)*180/math.pi
  if verbose:
    print "Point", pt
    print "pt bearing", pt.bearing
    print "view bearing", yaw
    print "diff", diff
    print
  return diff < 85

def picknearestll(dict2d, tag):
  """Returns closest value in dict"""
  locs = []
  for pixel, loc in dict2d.items():
    if loc:
      locs.append((tag.xydistance(loc), loc))
  locs.sort()
  return locs[0][1]

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

def distance3dt(d1, d2):
  """Finds distance between (lat,lon,alt) points."""
  return distance3d6(d1[0], d1[1], d1[2], d2[0], d2[1], d2[2])

def distance3d(d1, d2):
  """Finds distance between {'lat': deg, 'lon': deg, 'alt': meters} points."""
  x1, y1, z1 = d1['lat'], d1['lon'], d1['alt']
  x2, y2, z2 = d2['lat'], d2['lon'], d2['alt']
  return distance3d6(x1, y1, z1, z2, y2, z2)

def cameramat(wx, wy, fov):
    f = (wx-1) / (2*np.tan(fov/2))
    return np.array( [[ f, 0, (wx-1)/2. ],
                      [ 0, f, (wy-1)/2. ],
                      [ 0, 0, 1 ]]  )

def RfromYPR(yaw,pitch,roll):
    # relative yaw->pitch->roll
    yaw, pitch, roll = np.pi/180*yaw, np.pi/180*pitch, np.pi/180*roll
    Ry = np.array([[np.cos(yaw),0,np.sin(yaw)],
                   [0,1,0],
                   [-np.sin(yaw),0,np.cos(yaw)]])
    Rx = np.array([[1,0,0],
                   [0,np.cos(pitch),-np.sin(pitch)],
                   [0,np.sin(pitch),np.cos(pitch)]])
    Rz = np.array([[np.cos(roll),-np.sin(roll),0],
                   [np.sin(roll),np.cos(roll),0],
                   [0,0,1]])
    return np.dot(Ry,np.dot(Rx,Rz))

def YPRfromR(R):
    yaw   = 180 / np.pi * np.arctan2(R[0,2],R[2,2])
    pitch = 180 / np.pi * np.arcsin(-R[1,2])
    roll  = 180 / np.pi * np.arctan2(R[1,0],R[1,1])
    return yaw, pitch, roll

def xprodmat(x):
    return np.array([[0,-x[2],x[1]],
                     [x[2],0,-x[0]],
                     [-x[1],x[0],0]])

def NEtoLL(lat1, lon1, north, east):
    # computes lat2, lon2 which is a given distance east and north of lat1, lon1
    lat1 = math.pi/180 * lat1
    lon1 = math.pi/180 * lon1
    b = math.atan2(east,north)
    d = (north**2 + east**2)**0.5
    R = 6.3781e6
    lat2 = 180/math.pi * math.asin(math.sin(lat1)*math.cos(d/R) + math.cos(lat1)*math.sin(d/R)*math.cos(b))
    lon2 = 180/math.pi * ( lon1 + math.atan2(math.sin(b)*math.sin(d/R)*math.cos(lat1), math.cos(d/R)-math.sin(lat1)*math.sin(lat2)) )
    return lat2, lon2

def decomposeE(E):
    U,s,Vt = alg.svd(E)
    print 'singular values'
    print s
    S = np.array([[s[0],0,0],
                  [0,s[1],0],
                  [0,0,0]])
    W = np.array([[0,-1,0],
                  [1,0,0],
                  [0,0,1]])
    Z = np.array([[0,-1,0],
                  [1,0,0],
                  [0,0,0]])
    R = np.transpose(np.dot(U,np.dot(W,Vt))) # orientation of I2 in I1 frame
    y,p,r = YPRfromR(R)
    if abs(y)>90 or abs(p)>90 or abs(r)>90: # wrong "sign" on orientation
        R = np.transpose(np.dot(U,np.dot(np.transpose(W),Vt))) # orientation of I2 in I1 frame
    tx = np.dot(np.transpose(Vt),np.dot(Z,Vt))
    t = np.array([tx[2,1],tx[0,2],tx[1,0]])
    return R,t,s[0]/s[1]

def decomposeH(H):
    U,s,Vt = alg.svd(H)
    H = H / s[1]
    s = s / s[1]
    s2 = s**2
    print 'singular values'
    print s

    # decompose homography
    v1 = Vt[0,:]
    v2 = Vt[1,:]
    v3 = Vt[2,:]
    u1 = ( np.sqrt(1-s2[2])*v1 + np.sqrt(s2[0]-1)*v3 ) / np.sqrt(s2[0]-s2[2])
    u2 = ( np.sqrt(1-s2[2])*v1 - np.sqrt(s2[0]-1)*v3 ) / np.sqrt(s2[0]-s2[2])
    Hv2 = np.dot(H,v2)
    Hu1 = np.dot(H,u1)
    Hu2 = np.dot(H,u2)
    V2 = np.array( [ [ 0,-v2[2],v2[1] ],
                     [ v2[2],0,-v2[0] ],
                     [ -v2[1],v2[0],0 ] ] )
    HV2 = np.array( [ [ 0,-Hv2[2],Hv2[1] ],
                      [ Hv2[2],0,-Hv2[0] ],
                      [ -Hv2[1],Hv2[0],0 ] ] )
    U1 = np.transpose( np.array( [ v2,u1,np.dot(V2,u1) ] ) )
    U2 = np.transpose( np.array( [ v2,u2,np.dot(V2,u2) ] ) )
    W1 = np.transpose( np.array( [ Hv2,Hu1,np.dot(HV2,Hu1) ] ) )
    W2 = np.transpose( np.array( [ Hv2,Hu2,np.dot(HV2,Hu2) ] ) )

    # compute R and t
    R1 = np.dot(W1,np.transpose(U1))
    n1 = np.dot(V2,u1)
    n1 = n1 / alg.norm(n1)
    t1 = np.dot(R1,np.dot(H-R1,n1))
    t1 = t1 / alg.norm(t1) # position of camera 2 in camera 1 frame
    R1 = np.transpose(R1) # orientation of camera 2 in camera 1 frame
    y1,p1,r1 = YPRfromR(R1)
    print 'YPR, T, and N from homography: Candidate 1'
    print np.array([y1,p1,r1])
    print t1
    print n1

    R2 = np.dot(W2,np.transpose(U2))
    n2 = np.dot(V2,u2)
    n2 = n2 / alg.norm(n2)
    t2 = np.dot(R2,np.dot(H-R2,n2))
    t2 = t2 / alg.norm(t2) # position of camera 2 in camera 1 frame
    R2 = np.transpose(R2) # orientation of camera 2 in camera 1 frame
    y2,p2,r2 = YPRfromR(R2)
    print 'YPR, T, and N from homography: Candidate 2'
    print np.array([y2,p2,r2])
    print t2
    print n2

#    if abs(y1)>np.pi/2 or abs(p1)>np.pi/2 or abs(r1)>np.pi/2: # wrong orientation
    return R2,t2,n2
#    else:
#        return R1,t1,n1

# vim: et sw=2
