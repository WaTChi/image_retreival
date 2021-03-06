# as far as I can tell earthmine .info file related util

import os
import shutil
import math
import random
SIFTREGEXSTR = r'.*sift.txt$'
CHOGREGEXSTR = r'.*chog.txt$'
SURFREGEXSTR = r'.*surf.npy$'
IMGREGEXSTR = r'.*.jpg$'

def add_error((lat, lon), error):
    if error is 0:
        return (lat, lon)
    if type(error) is dict:
        random.seed(lat + lon + error['seed'])
        bearingDegrees = random.uniform(0, 360)
        error = abs(getattr(random, error['dist'])(*error['args']))**error['exponent']
    else:
        random.seed(lat + lon + 1)
        bearingDegrees = random.uniform(0, 360)
    print ">>> RANDOM BEARING >>> %d" % bearingDegrees
    dest = moveLocation(lat, lon, error, bearingDegrees)
    print (lat, lon), dest
    return dest

def moveLocation(lat1, lon1, d, bearingDegrees):
    """Returns a ddLocation that is d meters from view1
       along a great circle heading along the bearing"""

    R = 6.3781e6 #Earth's radius in meters
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)

    bearing = math.radians(bearingDegrees)
    lat2 = math.asin(math.sin(lat1)*math.cos(d/R) \
                     + math.cos(lat1)*math.sin(d/R)*math.cos(bearing))
    lon2 = lon1 + \
           math.atan2(math.sin(bearing)*math.sin(d/R)*math.cos(lat1),
                      math.cos(d/R) - math.sin(lat1)*math.sin(lat2))
    return math.degrees(lat2), math.degrees(lon2)

def getInfoCoord(path):
    """gets coordinates from .info file. makes assumption about format of file"""
    if os.path.exists(path):
        f = open(path, 'r')
        a=f.readline()
        lat=a.split("u'lat\': ")[1].split(",")[0]
        lon=a.split("u'lon\': ")[1].split("}")[0]
        return float(lat), float(lon)
    else:
        raise OSError("{p} does not exist.".format(p=path))
def getQueryCoord(fname):
    lat=fname.split(",")[1]
    lon=fname.split(",")[2][0:-4]
    return float(lat), float(lon)
def getQuerySIFTCoord(fname):
    #print fname.split(",")
    lat=fname.split(",")[0]
    lon=fname.split(",")[1][0:-13]
    return float(lat), float(lon)
getQueryCHOGCoord = getQuerySIFTCoord
getQuerySURFCoord = getQuerySIFTCoord
def getCellCoord(dname):
    lat, lon = dname.split(',')
    lat = float(lat)
    lon = float(lon)
    return float(lat), float(lon)
def getImgCoord(fname):
    """gets coordinates based on JPG file name. makes assumption about format of filename"""
    lat=fname.split(",")[0]
    lon=fname.split(",")[1][0:-9]
    return float(lat), float(lon)
def getSIFTCoord(fname):
    """gets coordinates based on SIFT file name. makes assumption about format of filename"""
    lat=fname.split(",")[0]
    lon=fname.split(",")[1][0:-13]
    return float(lat), float(lon)
getCHOGCoord = getSIFTCoord
getSURFCoord = getSIFTCoord
def getSIFTAngle(fname):
    """gets angle based on JPG file name. makes assumption about format of filename"""
    return int(fname[-12:-8])
def siftdistance(a, b):
  lat1, lon1 = getSIFTCoord(a)
  lat2, lon2 = getSIFTCoord(b)
  return distance(lat1, lon1, lat2, lon2)

def distance(lat1, lon1, lat2, lon2):
    """Gets geographic distance between locations using spherical law of cosines"""
    if lat1 == lat2 and lon1 == lon2:
        return 0
    R = 6.3781e6
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)
    p = min(1.0, math.sin(lat1)*math.sin(lat2) \
                  + math.cos(lat1)*math.cos(lat2) \
                  * math.cos(lon2 - lon1))
    d = math.acos(p)*R
    return d

def getbearing(lat1,long1, lat2, long2):
    """Takes two locations and returns the initial bearing (in degrees) from
       points1 to point2 where 0 is north."""
    dLon = math.radians(long2-long1)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)
    y = math.sin(dLon)*math.cos(lat2)
    x = math.cos(lat1)*math.sin(lat2) - \
        math.sin(lat1)*math.cos(lat2)*math.cos(dLon)
    return math.degrees(math.atan2(y,x))
    
def moveCellsIf(path, lat, lon, radius):
    if not os.path.exists(path):
        return
    dirs = os.listdir(path)
    dirs = [ dir for dir in dirs if os.path.isdir(os.path.join(path,dir))]
    for dir in dirs:
        # if l>100:
        lat2, lon2 = dir.split(',')
        lat2 = float(lat2)
        lon2 = float(lon2)
        if distance(lat,lon,lat2,lon2) < radius:
            dirpath = os.path.join(path, dir)
            outpath = os.path.join(path, lat,',',lon,',',radius)
            shutil.move(dirpath, outpath)
