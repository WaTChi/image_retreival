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

def get_spots_in_square(lat, lon, halflength, x, y):
    #fetches spots within square centered at lat,lon and halflength given.
    center = earthMine.ddLocation(lat, lon)    
    print "Getting list of views within {rad} meters...".format(rad=halflength),
    conn = earthMine.ddObject()
    r=100
    if halflength <= r:
        views = earthMine.getFrontalViews(conn, lat, lon, halflength, maxResults=1000, width=x, height=y)
    else:
        #Generate other centers...
        #Throw down a square, sprinkle centers over it...
        maxLat = earthMine.moveLocation(center, halflength, 0)["lat"]
        minLon = earthMine.moveLocation(center, halflength, 270)["lon"]

        N = int(math.ceil(halflength/r))
        M = N-1
        centers=[]
        #Major Grid = upper left corner
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

        #Process list, kill duplicates
        #Easiest way: two lists, one without dupes.
        views = []
        for v in iviews:
            #Not a duplicate?
            dup=False
            for i in views:
                if v["view-location"]== i["view-location"]:
                    dup=True
                    break
            if not dup:
                views.append(v)
    print "Found {0} bubbles...".format(len(views))
    return views

    
def get_cell(lat, lon, halflength, outdir, imagesPerBubble=12):
    
    spots = get_spots_in_square(lat, lon, halflength, x=earthMine.DEFWIDTH, y=earthMine.DEFHEIGHT)
    #Generate directory from LAT, LON, r
    DIRNAME = os.path.join(outdir, str(lat)+"," + str(lon)+"," + str(halflength))
    print "Making directory: "+DIRNAME
    if not os.path.exists(DIRNAME):
        try:
            os.makedirs(DIRNAME)
        except Exception:
            print "Error making directory...quitting..."
            return
    print "Made directory: "+DIRNAME
    #Write out all the view info so we can go get things later if need be:
    f = open(os.path.join(DIRNAME,"spots.bin"), "wb")
    cPickle.dump(spots, f)
    f.close()
    print "Wrote {0} views to {1}...".format(len(spots), os.path.join(DIRNAME,"views.bin"))
    conn = earthMine.ddObject()
    for spot in spots:
        views = earthMine.buildCylinderFromView(conn, spot, imagesPerBubble)
        earthMine.saveViews(views, DIRNAME)
    print "Images are stored in ./"+DIRNAME
    return
    

halflength = 100
x = earthMine.DEFWIDTH
y = earthMine.DEFHEIGHT
getDepth = False

if __name__=="__main__":

    if len(sys.argv) < 3:
        print "USAGE: {prog} lat lon halflength outdir]".format(prog=sys.argv[0])
        sys.exit()
    

#    #SYSTEM PARAMETERS
    imagesPerCylinder = 12
    lat = float(sys.argv[1])
    lon = float(sys.argv[2])
    halflength = float(sys.argv[3])
    outdir = sys.argv[4]
    get_cell(lat, lon, halflength, outdir, imagesPerCylinder)