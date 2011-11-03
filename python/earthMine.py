# An implementation of the EarthMine DirectData API using their
# HTTP Request Protocol

import httplib2
import urllib
import time
import json
import hashlib
import math
import numpy
import os
from PIL import Image

#DD_KEY="v12ypuyn9a2renzw9ivvu5fg"
#DD_SECRET="Q4xfUF0vqO"
#Joe's
DD_KEY = "bv5tdq3uw9mjudeqf6kjctm3"
DD_SECRET = "M3jXvKKTKg"
DD_MAX_IMAGE_SIZE = 2048
# DEFWIDTH = 512
# DEFHEIGHT = 512
#DEFWIDTH = 768
#DEFHEIGHT = 512
DEFWIDTH = 2048
DEFHEIGHT = 1371
DEFFOV = 60
DD_DELAY = 0.20
DD_VIEWPIXELMAX = 400
DD_VIEWPIXELSIDE = 20
EARTHMINE_URL = "http://cloud.earthmine.com/service"

DEBUG=0

class ddError(Exception):
    """Default Exception class for ddErrors."""
    def __init__(self, status, reason):
        self.status = status
        self.reason = reason

    def __str__(self):
        return "STATUS: "+repr(self.status)+" REASON:"+self.reason

class ddLocation( dict ):
    def __init__(self, lat, lon, alt=None):
        dict.__init__(self)
        #check safety
        try:
            if lat > 90.0 or lat < -90.0:
                raise ddError(1, "Latitude out of range [90.0,-90.0]")
            elif lon > 180.0 or lon < -180.0:
                raise ddError(2, "Longitude out of range [180.0,-180.0]")
            elif alt is not None and alt < 0:
                raise ddError(3, "Altitude less than zero")
        except ValueError:
            raise ddError(4, "ddLocation initialized with non-numerical parameter")
        else:
            self["lat"] = lat
            self["lon"] = lon
            if alt is not None:
                self["alt"] = alt

    def __setitem__(self, key, value):
        if key not in ["lat", "lon", "alt"]:
            raise ddError(5, "Can only modify lat, lon, alt")
        elif key is "lat" and (value > 90.0 or value < -90.0):
            raise ddError(1, "Latitude out of range [90.0,-90.0]")
        elif key is "lon" and (value > 180.0 or value < -180.0):
            raise ddError(2, "Longitude out of range [180.0,-180.0]")
        elif key is "alt" and value < 0:
            raise ddError(3, "Altitude less than zero")
        else:
            dict.__setitem__(self, key, value)

class ddImageSize( dict ):
    def __init__(self, width, height):
        dict.__init__(self)
        try:
            if width < 0 or width > DD_MAX_IMAGE_SIZE:
                raise ddError(5, "Image width out of bounds [0,{max}]".format(max=DD_MAX_IMAGE_SIZE))
            elif height < 0 or height > DD_MAX_IMAGE_SIZE:
                raise ddError(6, "Image height out of bounds [0,{max}]".format(max=DD_MAX_IMAGE_SIZE))
        except ValueError:
            raise ddError(7, "ddImageSize initialized with bad parameter")
        else:
            self["width"]=width
            self["height"]=height

    def __setitem__(self, key, value):
        if key not in ["width", "height"]:
            raise ddError(5, "Can only modify width, height")
        elif key is "width" and (value < 0 or value > DD_MAX_IMAGE_SIZE):
            raise ddError(5, "Image width out of bounds [0,{max}]".format(max=DD_MAX_IMAGE_SIZE))
        elif key is "height" and (value < 0 or value > DD_MAX_IMAGE_SIZE):
            raise ddError(6, "Image height out of bounds [0,{max}]".format(max=DD_MAX_IMAGE_SIZE))
        else:
            dict.__setitem__(self, key, value)

class ddView( dict ):
    def __init__(self, imageSize, fieldOfView, isKnownOccluded,
                 url, viewLocation, viewId):

        dict.__init__(self)
        self["image-size"] = imageSize
        self["view-location"] = viewLocation
        self["field-of-view"] = fieldOfView
        self["is-known-occluded"] = isKnownOccluded
        self["url"] = url
        self["id"] = viewId

    def __setitem__(self, key, value):
        if key not in ["image-size", "view-location", "field-of-view",
                       "url", "view-location", "id", "is-known-occluded"]:
            raise ddError(0, "Fail")
        else:
            dict.__setitem__(self, key, value)

class ddViewRequest( dict ):
    def __init__(self, imageSize, viewSubject, FOV,
                  maxResults = 12, searchRadius = 60):
        dict.__init__(self)
        self["image-size"] = imageSize
        self["view-subject"]= viewSubject
        self["field-of-view"] = FOV
        self["max-results"]=maxResults
        self["search-radius"]=searchRadius

    def __setitem__(self, key, value):
        if key not in ["image-size", "view-subject", "field-of-view",
                        "max-results", "search-radius"]:
            raise ddError(0, "Fail")
        else:
            dict.__setitem__(self, key, value)

class ddViewPixel( dict ):
    def __init__(self, x, y):
        dict.__init__(self)
        self["x"] = x
        self["y"] = y

    def __setitem__(self, key, value):
        if key not in ["x", "y"]:
            raise ddError(0, "Fail")
        else:
            dict.__setitem__(self, key, value)

    def imprt( data ):
        """Imports a dictionary with the right structure
           as a ddViewPixel object."""
        try:
            self["x"] = data["x"]
            self["y"] = data["y"]
        except KeyError:
            raise ddError(32, "Import failed because of incorrect dictionary.")

class ddViewLocationForViewPixelsRequest( dict ):
    def __init__(self, viewID, viewPixels):
        """viewID should be the earthmine view-id string.
           viewPixels should be a list of viewPixel objects."""
        dict.__init__(self)
        self["view-id"] = viewID
        self["view-pixels"] = viewPixels

    def __setitem__(self, key, value):
        if key not in ["view-id", "view-pixels"]:
            raise ddError(0, "Fail")
        else:
            dict.__setitem__(self, key, value)

class ddObject():
    """Handles connections to the earthmine Direct Data API.
       Supports all Direct Data Request types.

       General workflow:
          ddObject.buildXYZRequest(PARAMS)
          try:
              content = ddObject.sendRequest()
          except ddError:
              HANDLE"""

    def __init__(self):
        #Construct request headers
        self.requestURL = EARTHMINE_URL
        self.httpObject = httplib2.Http()
        self.headers={"x-earthmine-auth-id": DD_KEY}
        self.JSONData = ""
        self.JSONEncoder = json.JSONEncoder(ensure_ascii=False)
        self.JSONDecoder = json.JSONDecoder()

    def signRequest(self):
        currTime = int(time.time())
        sig = hashlib.md5(DD_KEY+DD_SECRET+str(currTime)).hexdigest()
        #fullURL = self.requestURL + "?sig="+sig
        fullURL = self.requestURL + "?sig="+sig+"&timestamp="+str(currTime)
#        print fullURL
        return fullURL

    def buildGetViews(self, width, height, FOV=DEFFOV,
                         lat=0.0, lon=0.0, alt=None,
                         maxResults=12, searchRadius=60.0):
        requestDict={"operation": "get-views",
                     "parameters": {
                         "request": ddViewRequest(ddImageSize(width, height),
                                                  ddLocation(lat, lon, alt),
                                                  FOV, maxResults, searchRadius)

                         }
                     }
        self.JSONData = self.JSONEncoder.encode(requestDict)

    def buildGetPano(self, lat=0.0, lon=0.0, maxResults=12, searchRadius=60.0):
        requestDict={"operation": "get-panoramas",
                        "parameters": {
                            "request":{
                                "subject-location":{
                                    "lat": lat,
                                    "lon": lon
                                    },
                             "max-results": maxResults,
                             "search-radius": searchRadius
                             }
                        }
                    }
        self.JSONData = self.JSONEncoder.encode(requestDict)

    def processViewSearchResult():
        pass

    def adjustView(self, view, FOV=DEFFOV, pan=0, tilt=0, width=False, height=False):
        """Adjust any view passed to the object and return the new view.
           The old view is not kept. Pan and tilt are in degrees"""
        if width == False:
            width = view["image-size"]["width"]
        if height == False:
            height = view["image-size"]["height"]
        request = {"operation":"adjust-view",
                   "parameters": {
                       "request": {
                           "view-id": view["id"],
                           "field-of-view":FOV,
                           "image-size":{"width":width, "height":height},
                           "pan":pan,
                           "tilt": tilt
                           }
                       }
                   }
        self.JSONData = self.JSONEncoder.encode(request)
        return self.sendRequest()

    def getdepthurl(self, view):
        """Retrieve a url of a depth image for the given view."""
        request = {"operation":"get-depth-image",
                   "parameters" : {
                       "request" : {
                           "view-id": view["id"]
                           }
                       }
                   }
        self.JSONData = self.JSONEncoder.encode(request)
        return self.sendRequest()

    def getLocationsForViewPixels( self, viewID, viewPixels):
        request = {
            "operation" : "get-locations-for-view-pixels",
            "parameters" : {
                "request" :
                    ddViewLocationForViewPixelsRequest(viewID, viewPixels)

                }
            }

        self.JSONData = self.JSONEncoder.encode(request)
        response = self.sendRequest()
        #correllate the pixels and locations in a dictionary
        return [(request["parameters"]["request"]["view-pixels"][i], response["locations"][i]) for i in range(len(response["locations"]))]

    def processResult(self, content):
        """Generic Handler for processing all earthmine repsonses.
           Checks for internal exceptions and returns a dictionary of
           the result."""

        resultDict = self.JSONDecoder.decode(content)
        #check for errors
        if resultDict["exception"] is None:
            return resultDict["result"]
        else:
            raise ddError(13,
                             "Exception in response: {resp}".format(resp=
                              repr(resultDict["exception"])))


    def sendRequest(self):
        if self.JSONData == "":
            pass
        else:
            fullURL = self.signRequest()
            self.httpObject = httplib2.Http()
            try:
                resp, content = self.httpObject.request(fullURL, 'POST',
                                                        headers=self.headers,
                                                        body=self.JSONData)
            except: #Try again for one error
                resp, content = self.httpObject.request(fullURL, 'POST',
                                                        headers=self.headers,
                                                        body=self.JSONData)
            #process header errors here
        while True:
            if resp.status == 200:
                return self.processResult(content)
            elif resp.status == 504:
                resp, content = self.httpObject.request(fullURL, 'POST',
                                                    headers=self.headers,
                                                    body=self.JSONData)
            else:
                raise ddError(resp.status, resp.reason + content)

#Higher Level functions
def LocationSubtract(view1, view2):
##    print "v1: ",
##    print view1
##    print "v2: ",
##    print view2
    """Gets distance between locations using spherical law
       of cosines"""
    R = 6371e3
    try:
        lat1 = math.radians(view1["lat"])
        lon1 = math.radians(view1["lon"])
    except KeyError:
        raise ddError(3, "view1 did not have members 'lat' or 'lon'")

    try:
        lat2 = math.radians(view2["lat"])
        lon2 = math.radians(view2["lon"])
    except KeyError:
        raise ddError(3, "view1 did not have members 'lat' or 'lon'")

    d = math.acos(math.sin(lat1)*math.sin(lat2) \
                  + math.cos(lat1)*math.cos(lat2) \
                  * math.cos(lon2 - lon1)) * R
    return d

#def LocationSubtractHaver(view1, view2):
#    """Subtracts two ddLocations and returns the distance between them
#       using the Haversine formula"""
#
#    R = 6371e3 #Earth's radius in meters
#    dLat = math.radians(view1["lat"]-view2["lat"])
#    dLon = math.radians(view1["lon"]-view2["lon"])
#    a = math.sin(dLat/2)*math.sin(dLat/2) \
#        + math.cos(math.radians(view1["lat"])) \
#        * math.cos(math.radians(view2["lat"])) \
#        * math.sin(dLon/2)*math.sin(dLon/2)
#    c = math.atan2(math.sqrt(a), math.sqrt(1-a))
#    return R * c

#def tobearing(rad):
#    return (math.degrees(rad) + 360) % 360

def getbearing(view1, view2):
    """Takes two ddLocations and returns the initial bearing (in degrees) from
       view1 to view2 where 0 is north."""
    dLon = math.radians(view2["lon"]-view1["lon"])
    lat1 = math.radians(view1["lat"])
    lat2 = math.radians(view2["lat"])
    y = math.sin(dLon)*math.cos(lat2)
    x = math.cos(lat1)*math.sin(lat2) - \
        math.sin(lat1)*math.cos(lat2)*math.cos(dLon)
    return math.degrees(math.atan2(y,x))
#    return tobearing(math.atan2(y,x))



def moveLocation(view1, d, bearingDegrees):
    """Returns a ddLocation that is d meters from view1
       along a great circle heading along the bearing"""

    R = 6371e3 #Earth's radius in meters
    try:
        lat1 = math.radians(view1["lat"])
        lon1 = math.radians(view1["lon"])
    except KeyError:
        raise ddError(3, "Passed view was not a ddLocation...")

    bearing = math.radians(bearingDegrees)
    lat2 = math.asin(math.sin(lat1)*math.cos(d/R) \
                     + math.cos(lat1)*math.sin(d/R)*math.cos(bearing))
    lon2 = lon1 + \
           math.atan2(math.sin(bearing)*math.sin(d/R)*math.cos(lat1),
                      math.cos(d/R) - math.sin(lat1)*math.sin(lat2))

    return ddLocation(math.degrees(lat2), math.degrees(lon2))

def moveLocation4(lat1, lon1, d, bearingDegrees):
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

def getCell(ddConnection, lat, lon, radius):
    conn = ddObject()
    views = earthMine.ddGetViews(conn, lat, lon, radius, maxResults=100)


def ddGetViews(ddConnection, lat, lon, radius=60, maxResults=12, FOV=DEFFOV, width=DEFWIDTH, height=DEFHEIGHT):
    """Gets a list of views within RADIUS meters of
       the given latitude and longitude."""

    try:
        ddConnection.buildGetViews(width, height, searchRadius=radius, maxResults=maxResults, FOV=FOV, lat=lat, lon=lon)
        if DEBUG > 0:
            print ddConnection.JSONData
        result = ddConnection.sendRequest()
    except Exception as msg:
        print "sent: {0}".format(ddConnection.JSONData)
        raise ddError(0, "ddConnection object raised an exception: {msg}".format(msg=msg))

    try:
        return result["views"]
    except KeyError:
        return None

def ddGetPanos(ddConnection, lat, lon, radius=60, maxResults=12):
    """Gets a list of panos within RADIUS meters of
       the given latitude and longitude."""
    try:
        ddConnection.buildGetPano(searchRadius=radius, maxResults=maxResults, lat=lat, lon=lon)
        if DEBUG > 0:
            print ddConnection.JSONData
        result = ddConnection.sendRequest()
    except Exception as msg:
        print "sent: {0}".format(ddConnection.JSONData)
        raise ddError(0, "ddConnection object raised an exception: {msg}".format(msg=msg))

    try:
        return result["panoramas"]
    except KeyError:
        return None

def getFrontalViews(ddConnection, lat, lon, radius=60, maxResults=12, FOV=DEFFOV, width=DEFWIDTH, height=DEFHEIGHT):
    panos = ddGetPanos(ddConnection, lat, lon, radius=radius, maxResults=maxResults)
    views = [getFrontalView(ddConnection, pano, FOV, width, height) for pano in panos]
    return views

def getFrontalView(ddConnection, pano, FOV, width, height):
    lat, lon = moveLocation4(pano['location']['lat'], pano['location']['lon'], .25, pano['pano-orientation']['yaw'])
    return ddGetViews(ddConnection, lat, lon, 10, maxResults=1, FOV=FOV, width=width, height=height)[0]

def saveViews(views, outputDir, getDepth=False):
    if not os.path.exists(outputDir):
            try:
                os.makedirs(outputDir)
            except Exception:
                print "Error making directory...quitting..."
                return
    count=0
    for view in views:
        error = 0
        while error < 10:
            try:
                fname = os.path.join(outputDir, str(view["view-location"]["lat"])+","+str(view["view-location"]["lon"])+"-"+"{0:04}".format(count))
                urllib.urlretrieve(view["url"]["href"], fname+".jpg")

                #Get depth image
                #Other depth things...
                if getDepth:
                    conn = earthMine.ddObject()
                    locs = earthMine.ddGetImageLocs(view, conn)
                    ThreeDData = earthMine.ddImageLocstoLPT(view, locs)
                    earthMine.ddWriteLocsFile(locs, fname+".locs")
                    earthMine.ddWriteLPTFile(ThreeDData, fname+".lpt")
                    earthMine.ddWriteDepthFile(ThreeDData, fname+".depth")
                    earthMine.ddMakeDepthImage(view, ThreeDData, fname+"depth.tif")
                error = 10
            except IOError:
                error = error + 1

        f = open(fname+".info", "w")
        f.write(repr(view))
        f.close()
        count += 1

def buildCylinderFromViewNoStreet(conn, spot, imagesPerCylinder):
    views = []
    FOV = max(360 / imagesPerCylinder + 20, 60) #For 4 views, FOV is 110, max allowed
    panAmount = 360 / imagesPerCylinder / 2
    for i in range(1, imagesPerCylinder):
        if(i != imagesPerCylinder/2):
            time.sleep(DD_DELAY)
            views.append(conn.adjustView(spot, FOV, panAmount*i)["view"])
    return views

def buildCylinderFromView(conn, spot, imagesPerCylinder):
    views = []
    FOV = max(360 / imagesPerCylinder + 20, 60) #For 4 views, FOV is 110, max allowed
    panAmount = 360 / imagesPerCylinder / 2
    # conn = earthMine.ddObject()
    # count = 0
    # for spot in spots:
        # dbgmsg("Cylinder {n} of {k}".format(n=count, k=len(spots)))
    views.append(spot)
    for i in range(1, imagesPerCylinder):
        time.sleep(DD_DELAY)
        views.append(conn.adjustView(spot, FOV, panAmount*i)["view"])
        # count += 1
    return views

def ddMakeImageCyl(ddConnection, lat, lon, numImages, width=DEFWIDTH, height=DEFHEIGHT):
    """Given a latitude and longitude, get closest view
       and derive an image polygon from around this view. FOV
       is set automaticallty to overlap views in the polygon.
       The number of images must be at least 4"""

    if numImages+1-1 != numImages:
        raise ddError(0, "numImages must be a number")
    if numImages < 4:
        raise ddError(0, "Not enough views specified")
    views = []
    FOV = max(360 / numImages + 20, 60) #For 4 views, FOV is 110, max allowed
    panAmount = 360 / numImages / 2 #Some unkown bug makes this division necessary
    nView = ddGetViews(ddConnection, lat, lon, FOV=FOV, width=width, height=height)[0]
    views.append(nView)

    for i in range(1, numImages):
        time.sleep(DD_DELAY) #Don't exceed my rate limit
        views.append(ddConnection.adjustView(views[i-1], FOV, panAmount)["view"])

    return views

class coord(tuple):
    def __new__(self,x=0,y=0,z=0):
        self = (x,y,z)


def ddGetImageLocs(view, conn):
    """Given an earthmine view, get 3D points for every pixel in the image.
       the camera is assumed to sit at world point (0,0,0)."""
    x = 0
    y = 0
    locs = []
    while x < view["image-size"]["width"]:
        while y < view["image-size"]["height"]:
            xmax = min(x+DD_VIEWPIXELSIDE, view["image-size"]["width"])
            ymax = min(y+DD_VIEWPIXELSIDE, view["image-size"]["height"])
            pixels = [ddViewPixel(i,j) for i in range(x,xmax)
                      for j in range(y,ymax)]
            print("x={0}, y={1}, {2} pixels.".format(x, y, len(pixels)))
            locs.extend(conn.getLocationsForViewPixels(view["id"],pixels))
            y = y+DD_VIEWPIXELSIDE
        x = x+DD_VIEWPIXELSIDE
        y = 0
    return locs

def ddImageLocstoLPT(view, locs):
    # Translate location for each pixel into x,y,z where (0,0,0) is the camera center.
    # First, convert ECEF-geodetic to ECEF-rectangular
    originR = ddViewLocationToECEFR(view["view-location"])
    #Compute rotations and translations to center origin at 0 and make a
    #tangent plane aligned with east and north.
    #A good reference can be found at:
    #http://psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf
    t = originR
    theta = math.radians(view["view-location"]["lat"])
    phi = math.radians(view["view-location"]["lon"])
    R = numpy.matrix([[-1*math.sin(phi), math.cos(phi), 0],
                     [-1*math.cos(phi)*math.sin(theta), -1*math.sin(phi)*math.sin(theta), math.cos(theta)],
                     [math.cos(phi)*math.cos(theta), math.cos(theta)*math.sin(phi), math.sin(theta)]])
    ThreeDData = []
    for loc in locs: #Process image points into 3D points.
        if loc[1] is not None:
            ECEFPt = ddViewLocationToECEFR(loc[1])
            LPTPt = R * (ECEFPt - t)
            ThreeDData.append( ((loc[0]["x"],loc[0]["y"]),
                           LPTPt))
        else:
            ThreeDData.append( ((loc[0]["x"],loc[0]["y"]), None))


    return ThreeDData

def ddMakeDepthImage(view, ThreeDData, filename):
    im = Image.new("I", (view["image-size"]["width"],
                         view["image-size"]["height"]))
    for p in ThreeDData:
        if p[1] is not None:
            a = p[1]
            pix = numpy.sqrt(numpy.dot(a.T,a))*100 #cm
            im.putpixel(p[0], pix)
	else:
            im.putpixel(p[0], 0)
    im.save(filename)

def ddWriteLocsFile(locs, filename):
    f = open(filename,"w")
    for l in locs:
        if l[1] is not None:
            f.write("{0}, {1}, {2}, {3}, {4}\n".format(l[0]["x"],l[0]["y"], l[1]["lat"], l[1]["lon"], l[1]["alt"]))
        else:
            f.write("{0}, {1}, NULL, NULL, NULL\n".format(l[0]["x"], l[0]["y"]))
    f.close()

def ddWriteLPTFile(ThreeDData, filename):
    f = open(filename,"w")
    for d in ThreeDData:
        if d[1] is not None:
            f.write("{0}, {1}, {2}, {3}, {4}\n".format(d[0][0], d[0][1], d[1][0], d[1][1], d[1][2]).replace("[[","").replace("]]",""))
        #else:
            #f.write("{0}, {1}, {2}\n".format(d[0][0], d[0][1], "NULL"))
    f.close()

def ddWriteDepthFile(ThreeDData, filename):
    f = open(filename,"w")
    for d in ThreeDData:
        if d[1] is not None:
            f.write("{0}, {1}, {2}\n".format(d[0][0], d[0][1], numpy.sqrt(numpy.dot(d[1].T,d[1]))*100).replace("[[","").replace("]]",""))
        else:
            f.write("{0}, {1}, 0\n".format(d[0][0], d[0][1]))
    f.close()

def ddViewLocationToECEFR(viewLocation):
    """returns earth-centered, earth-fixed, rectangular coordinates"""
    a = 6378137.0 #Equatorial Radius
    b = 6356752.3 #Polar Radius
    theta = math.radians(viewLocation["lat"])
    R = math.sqrt( ((a*a*math.cos(theta))**2 + (b*b*math.sin(theta))**2) /
                   ((a*math.cos(theta))**2 + (b*math.sin(theta))**2) )
    r = R+viewLocation["alt"] #Origin R in meters.
    phi = math.radians(viewLocation["lon"])
    return numpy.matrix([r*math.cos(phi)*math.cos(theta),
            r*math.sin(phi)*math.cos(theta),
            r*math.sin(theta)]).T

def ddMakeDepthMap(lat, lon, conn, name, width=DEFWIDTH, height=DEFHEIGHT):
    """Given an earthmine view, build a sparse
    point cloud of depth around the center of the 'bubble'."""

    #Make a square of images that encircle the current view.
    views = ddMakeImageCyl(conn, lat, lon, 12, width, height)

    print "Done getting views!"
    #Get depth map for each view
    count = 0
    for view in views:
        locs = ddGetImageLocs(view, conn)
        ThreeDData = ddImageLocstoLPT(view, locs)
        ddWriteLocsFile(locs, name+repr(count)+".locs")
        ddWriteLPTFile(ThreeDData, name+repr(count)+".lpt")
        ddWriteDepthFile(ThreeDData, name+repr(count)+".depth")
        ddMakeDepthImage(view, ThreeDData, name+repr(count)+"depth.tif")
        urllib.urlretrieve(view["url"]["href"], name+repr(count)+".jpg")
        count = count + 1
    return

def ddGetAllPixels(pixels, viewId, keep_None=False):
    """fetches an arbitrary amount of pixels from EarthMine Direct Data"""
    conn = ddObject()
    viewPixels = [ddViewPixel(p[0], p[1]) for p in pixels]
    locs = {}
    while viewPixels:
        response = None
        for retry in range(3):
            try:
                if retry:
                    print 'try %d' % retry
                response = conn.getLocationsForViewPixels(viewId, viewPixels[:490])
                break
            except Exception, e:
                print e
        if response is None:
            raise Exception
        viewPixels = viewPixels[490:] # limit for api
        for pixel, loc in response:
            if loc or keep_None: # has valid mapping
                locs[(pixel['x'], pixel['y'])] = loc
    return locs # map of (x,y) to coords
