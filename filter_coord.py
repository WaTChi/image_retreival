import math
import re
import sys
import os
import shutil

def filter(resultFile, imgDir):
    """moves match results specified in resultFile to a match directory from specified imgDir"""
    try:
        files = os.listdir(imgDir)
    except OSError:
        raise
    if os.path.exists(resultFile):
        cellnum=resultFile.split('.match')[0].split('results')[1]
        f = open(resultFile, 'r')
        for line in f:
            a = line.split()
            queryfile = a[0]
            matches = a[1:len(a)]
            odir = queryfile.split('sift.txt')[0] + '-matches_in_cell'+cellnum
            if not os.path.exists(odir):
                os.mkdir(odir)
            for match in matches:
                fname = "image%04d" % int(match)
                results = [f for f in files if re.search(fname, f)]
                for result in results:
                    shutil.copy(imgDir+"\\"+result, odir)
                # [lat, lon] = info.getCoord(cellDir+fname)

if __name__=="__main__":
    print len(sys.argv)
    if len(sys.argv) != 3:
        print "USAGE: {0} MatchFile cellDir".format(sys.argv[0])
        sys.exit
    else:
        filter(sys.argv[1], sys.argv[2])