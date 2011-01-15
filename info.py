import os
import shutil
import math
SIFTREGEXSTR = r'.*sift.txt$'
CHOGREGEXSTR = r'.*chog.txt$'
IMGREGEXSTR = r'.*.jpg$'
  
def moveLocation(lat1, lon1, d, bearingDegrees):
    """Returns a ddLocation that is d meters from view1
       along a great circle heading along the bearing"""

    R = 6371e3 #Earth's radius in meters
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
    lat=fname.split(",")[1]
    lon=fname.split(",")[2][0:-8]
    return float(lat), float(lon)
getQueryCHOGCoord = getQuerySIFTCoord
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
def getSIFTAngle(fname):
    """gets angle based on JPG file name. makes assumption about format of filename"""
    return int(fname[-12:-8])
def siftdistance(a, b):
  lat1, lon1 = getSIFTCoord(a)
  lat2, lon2 = getSIFTCoord(b)
  return distance(lat1, lon1, lat2, lon2)

def distance(lat1, lon1, lat2, lon2):
    """Gets distance between locations using spherical law of cosines"""
    if lat1 == lat2 and lon1 == lon2:
        return 0
    R = 6371e3
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)
    d = math.acos(math.sin(lat1)*math.sin(lat2) \
                  + math.cos(lat1)*math.cos(lat2) \
                  * math.cos(lon2 - lon1)) * R
    return d
    
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
