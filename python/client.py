import info
import util
import os
import corr
import query
import time
from querySystemCopy import combine_ransac, Img, draw_top_corr
from config import INFO

params = query.PARAMS_DEFAULT.copy()
params.update({
  'checks': 1024,
  'trees': 1,
  'vote_method': 'matchonce',
})
cellradius = 236.6
ambiguity = 50
matchdistance = 25
ncells = 7
maindir = os.path.expanduser('~/shiraz')
#maindir = os.path.expanduser('~etzeng/shiraz')
dbdir = os.path.join(maindir, 'Research/cells/g=100,r=d=236.6/')
#matchdir = os.path.expanduser('~/results/%s' % query.searchtype(params))
matchdir = "/tmp/results/%s" % query.searchtype(params)
dbdump = os.path.join(maindir, "Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829")
if not os.path.isdir(matchdir):
    os.makedirs(matchdir)

def match(siftfile, imagefile, lat=None, lon=None):
    querydir = os.path.dirname(siftfile)
    siftfile = os.path.basename(siftfile)
    if lat == None or lon == None:
        lat, lon = info.getQuerySIFTCoord(siftfile)
    closest_cells = util.getclosestcells(lat, lon, dbdir)
    outputFilePaths = []
    cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:ncells] if dist < cellradius + ambiguity+matchdistance]
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(lat, lon, latcell, loncell)
        outputFilePath = os.path.join(matchdir, siftfile + ',' + cell + ',' + str(actualdist)  + ".res")
        outputFilePaths.append(outputFilePath)
    # start query
    timer = time.time()
    query.run_parallel(dbdir, [c for c,d in cells_in_range], querydir, siftfile, outputFilePaths, params)
    INFO("--> Running query: " + str(time.time()-timer) + "s")
    # end query
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(lat, lon, latcell, loncell)
        outputFilePath = os.path.join(matchdir, siftfile + ',' + cell + ',' + str(actualdist)  + ".res")
    combtime = time.time()
    comb_matches = corr.combine_matches(outputFilePaths)
    INFO("--> COMB: " + str(time.time() - combtime) + "s")
    rtime = time.time()
    combined = combine_ransac(comb_matches)
    INFO("--> RANSAC: " + str(time.time() - rtime) + "s")
    topentry = combined[0]
    matchedimg = topentry[0]
    matches = comb_matches[matchedimg + 'sift.txt']
    return matchedimg, matches

def draw_corr(queryimgpath, matchedimg, matches, matchoutpath=None):
    F, inliers = corr.find_corr(matches)
    matchimgpath = os.path.join(dbdump, '%s.jpg' % matchedimg)
    if matchoutpath == None:
        matchoutpath = os.path.expanduser('~/client-out.png')
    corr.draw_matches(matches, queryimgpath, matchimgpath, matchoutpath, inliers)
    return F, inliers

def preprocess_image(inputfile, outputfile=None, width=768, height=512):
    """Use the convert utility to preprocess an image."""
    timer = time.time()
    if outputfile == None:
        outputfile = inputfile.rsplit(".",1)[0] + ".pgm"
    os.system("convert {0} -resize {2}x{3} {1}".format(inputfile, outputfile, width, height))
    INFO("--> Image conversion: " + str(time.time()-timer) + "s")
    return outputfile

SIFTEXEC = os.path.join(maindir, 'Research/app/siftDemoV4/sift')

def extract_features(inputfile, outputfile=None, siftexec=SIFTEXEC):
    """Call the sift utility to extract sift features."""
    timer = time.time()
    if outputfile == None:
        outputfile = inputfile.rsplit(".",1)[0] + "sift.txt"
    os.system("{0} <{1} >{2}".format(siftexec, inputfile, outputfile))
    INFO("--> Feature extraction: " + str(time.time()-timer) + "s")
    return outputfile

if __name__ == '__main__':
    sift = os.path.expanduser('~/shiraz/DSC_7638,37.87162,-122.27223sift.txt')
    image = os.path.expanduser('~/shiraz/DSC_7638,37.87162,-122.27223.JPG')
    matchedimg, matches = match(sift, image)
    draw_corr(image, matchedimg, matches)

