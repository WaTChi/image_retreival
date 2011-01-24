import sys
import os
import earthMine
#import urllib
#import re
import math
import cPickle
import info
#import util

DEBUG = 0

def dbgmsg(string):
    if DEBUG > 0:
        print string

def createCells(lat, lon, radius, cellRad, uncertaintyRad, imgDir):
    center = earthMine.ddLocation(lat, lon)
    conn = earthMine.ddObject()
    cellDist=(3**.5) * (cellRad-uncertaintyRad)
    #Generate other centers...
    #Throw down a square of the radius, sprinkle centers over it...
    Lat1 = earthMine.moveLocation(center, radius, 0)["lat"]
    Lat2 = earthMine.moveLocation(center, radius, 180)["lat"]
    Lon1 = earthMine.moveLocation(center, radius, 90)["lon"]
    Lon2 = earthMine.moveLocation(center, radius, 270)["lon"]
    maxLat = max(Lat1, Lat2)
    minLat = min(Lat1, Lat2)
    maxLon = max(Lon1, Lon2)
    minLon = min(Lon1, Lon2)

    N = int(math.ceil(2*radius/cellDist))
    N2 = int(math.ceil(2*radius/((3**.5)*cellDist)))
    M = N-1
    M2 = N2-1
    centers=[]
    #Major Grid
    majstart = earthMine.ddLocation(maxLat, minLon)
#    minstart = earthMine.moveLocation(majstart, cellDist*math.sqrt(2), 135)
    minstart = earthMine.moveLocation(majstart, cellDist, 150)
    for i in range(N2):
        for j in range(N):
            centers.append(earthMine.moveLocation(majstart, cellDist*j, 90)) #Move east
        majstart = earthMine.moveLocation(majstart, (3**.5)*cellDist, 180) #Move south

    #Minor Grid
    for i in range(M2):
        for j in range(M):
            centers.append(earthMine.moveLocation(minstart, cellDist*j, 90))
        minstart = earthMine.moveLocation(minstart, (3**.5)*cellDist, 180)
        
    for c in centers:
        p=os.path.join(imgDir, str(c["lat"])+','+str(c["lon"]))
        util.copySIFTInRange(c["lat"], c["lon"], imgDir, p, cellRad)
    util.removeEmptyDir(imgDir)

def naiveHundredGet(lat, lon, radius, x, y):
    center = earthMine.ddLocation(lat, lon)    
    print "Getting list of views within {rad} meters...".format(rad=radius),
    conn = earthMine.ddObject()
    if radius <= 100:
        views = earthMine.getFrontalViews(conn, lat, lon, radius, maxResults=1000, width=x, height=y)
#        views = earthMine.ddGetViews(conn, lat, lon, radius, maxResults=100, FOV=60.0, width=x, height=y)
    else:
        #Generate other centers...
        #Throw down a square of the radius, sprinkle centers over it...
        Lat1 = earthMine.moveLocation(center, radius, 0)["lat"]
        Lat2 = earthMine.moveLocation(center, radius, 180)["lat"]
        Lon1 = earthMine.moveLocation(center, radius, 90)["lon"]
        Lon2 = earthMine.moveLocation(center, radius, 270)["lon"]
        maxLat = max(Lat1, Lat2)
        minLat = min(Lat1, Lat2)
        maxLon = max(Lon1, Lon2)
        minLon = min(Lon1, Lon2)

        r=100
#        N = 1+int(math.ceil((radius*2-100)/200))
#        if ((radius*2)/100)%2==0:
#            M = N-1            
#        else:
#            M = N
        N = int(math.ceil(radius/r))
        M = N-1            
        centers=[]
        #Major Grid
        majstart = earthMine.ddLocation(maxLat, minLon)
        minstart = earthMine.moveLocation(majstart, r*math.sqrt(2), 135)
        for i in range(N):
            for j in range(N):
                centers.append(earthMine.moveLocation(majstart, r*2*j, 90)) #Move east
            majstart = earthMine.moveLocation(majstart, r*2, 180) #Move south

        #Minor Grid
        for i in range(M):
            for j in range(M):
                centers.append(earthMine.moveLocation(minstart, r*2*j, 90))
            minstart = earthMine.moveLocation(minstart, r*2, 180)
            
        #DL all views
        iviews = []
        for c in centers:
            iviews.extend(earthMine.getFrontalViews(conn, c["lat"], c["lon"], r, maxResults=1000, FOV=60.0, width=x,height=y))

        #Process list, kill duplicates and views past the desired center...
        #Easiest way: two lists, one without dupes.
        views = []
        for v in iviews:
            #Inside bubble?
            if earthMine.LocationSubtract(center, v["view-location"]) < radius:
                #Not a duplicate?
                mindist=float("infinity")
                for i in views:
                    if v["view-location"]== i["view-location"]:
                        mindist=0
                if mindist!=0:
                    views.append(v)
#                dists = [earthMine.LocationSubtract(v["view-location"], i["view-location"]) for i in views]
#                if not dists or min(dists) > 4:
#                    views.append(v)
    print "Found {0} bubbles...".format(len(views))
    return views #ACTUAL END
    #return centers #DEBUG END
    
def getCell(lat, lon, radius, outdir, imagesPerBubble=12):
    
    spots = naiveHundredGet(lat, lon, radius, x=earthMine.DEFWIDTH, y=earthMine.DEFHEIGHT)
#    conn = earthMine.ddObject()
#    spots = earthMine.getFrontalViews(conn, lat, lon, radius, maxResults=1000, width=512, height=512)
    #Generate directory from LAT, LON
    DIRNAME = os.path.join(outdir, str(lat)+"," + str(lon))
    print "Making directory: "+DIRNAME
    if not os.path.exists(DIRNAME):
        try:
            os.makedirs(DIRNAME)
        except Exception:
            print "Error making directory...quitting..."
            return
    #Download images & depth info
    #Redo this to put each bubble in a separate directory!
    print "Made directory: "+DIRNAME
    #Write out all the view info so we can go get things later if need be:
    f = open(os.path.join(DIRNAME,"spots.bin"), "wb")
    cPickle.dump(spots, f)
    f.close()
    print "Wrote {0} views to {1}...".format(len(spots), os.path.join(DIRNAME,"views.bin"))
    conn = earthMine.ddObject()
    for spot in spots:
        views = earthMine.buildCylinderFromView(conn, spot, imagesPerBubble)
        #outdir = os.path.join(DIRNAME, str(views[0]["view-location"]["lat"]) + ", " + str(views[0]["view-location"]["lon"]))
        earthMine.saveViews(views, DIRNAME)
    print "Images are stored in ./"+DIRNAME
    return
    

radius = 100
x = earthMine.DEFWIDTH
y = earthMine.DEFHEIGHT
getDepth = False

if __name__=="__main__":

    if len(sys.argv) < 3:
        print "USAGE: {prog} lat lon radius outdir]".format(prog=sys.argv[0])
#        print "USAGE: {prog} lat lon [radius x y getDepth?]".format(prog=sys.argv[0])
        sys.exit()
    

#    #SYSTEM PARAMETERS
    imagesPerCylinder = 12
    lat = float(sys.argv[1])
    lon = float(sys.argv[2])
    radius = float(sys.argv[3])
    outdir = sys.argv[4]
    getCell(lat, lon, radius, outdir, imagesPerCylinder)
    

    #lattxt = sys.argv[1]
    #lontxt = sys.argv[2]
#    if len(sys.argv) >= 6:
#        x = int(sys.argv[4])
#        y = int(sys.argv[5])
#    if len(sys.argv) == 7:
#        getDepth = True
#
#    spots = naiveHundredGet(lat, lon, radius, x, y)
#
#
#    value = None
#    while value is None:
#        #value = raw_input("Getting {n} views...continue? (y/n)".format(
#            #n = len(views)))
#        value = "Y"
#        if value is "Y" or value is "y":
#            #Generate directory from LAT, LON
#            #DIRNAME = os.path.join("H:\\images\\", "img"+repr(lat)+"_"+repr(lon))
#            DIRNAME = os.path.join(OUTPATH, lattxt+","+lontxt)
#            print "Making directory: "+DIRNAME
#            try:
#                os.makedirs(DIRNAME)
#            except Exception:
#                print "Error making directory...quitting..."
#                sys.exit()
#            #Download images & depth info
#            #Redo this to put each bubble in a separate directory!
#            count = 1
#            print "Made directory: "+DIRNAME
#            #Write out all the view info so we can go get things later if need be:
#            f = open(os.path.join(DIRNAME,"spots.bin"), "wb")
#            cPickle.dump(spots, f)
#            f.close()
#            print "Wrote views to {0}...".format(os.path.join(DIRNAME,"views.bin"))
#            for y in spots:
#                conn = earthMine.ddObject()
#                views = earthMine.buildCylinderFromView(conn,y, imagesPerCylinder)
#                for x in views:
#                    print "Getting image {n}".format(n = count)
#                    error = 0
#                    while error < 10:
#                        try:
#                            fname = os.path.join(DIRNAME,"image"+"{0:04}".format(count))
#                            urllib.urlretrieve(x["url"]["href"], \
#                                       fname+".jpg")
#                                                   
#                            #Get depth image
#                            #Other depth things...
#                            if getDepth:
#                                conn = earthMine.ddObject()
#                                locs = earthMine.ddGetImageLocs(x, conn)
#                                ThreeDData = earthMine.ddImageLocstoLPT(x, locs)
#                                earthMine.ddWriteLocsFile(locs, fname+".locs")
#                                earthMine.ddWriteLPTFile(ThreeDData, fname+".lpt")
#                                earthMine.ddWriteDepthFile(ThreeDData, fname+".depth")
#                                earthMine.ddMakeDepthImage(x, ThreeDData, fname+"depth.tif")
#                            error = 10
#                        except IOError:
#                            error = error + 1
#                    
#                    f = open(fname+".info", "w")
#                    f.write(repr(x))
#                    f.close()
#                    count += 1
#
#            print "Images are stored in ./"+DIRNAME
#        elif value is "N" or value is "n":
#            print "Exiting..."
#            sys.exit()
#        else:
#            value = None
