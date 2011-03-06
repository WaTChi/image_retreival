import os
import corr
import info
import time
from config import *
import query
import queryContext as context

matchdir = '/tmp/client/%s' % query.searchtype(context.params)
if not os.path.exists(matchdir):
    os.makedirs(matchdir)
SIFTEXEC = os.path.join(context.maindir, 'Research/app/siftDemoV4/sift')
context.vars_init()

def preprocess_image(inputfile, outputfile=None, width=768, height=512):
    """Use the convert utility to preprocess an image."""
    timer = time.time()
    if outputfile == None:
        outputfile = inputfile.rsplit(".",1)[0] + ".pgm"
    os.system("convert {0} -resize {2}x{3} {1}".format(inputfile, outputfile, width, height))
    INFO("--> Image conversion: " + str(time.time()-timer) + "s")
    return outputfile

def extract_features(inputfile, outputfile=None, siftexec=SIFTEXEC):
    """Call the sift utility to extract sift features."""
    timer = time.time()
    if outputfile == None:
        outputfile = inputfile.rsplit(".",1)[0] + "sift.txt"
    os.system("{0} <{1} >{2}".format(siftexec, inputfile, outputfile))
    INFO("--> Feature extraction: " + str(time.time()-timer) + "s")
    return outputfile

def match(siftfile, lat=None, lon=None):
    if lat == None or lon == None:
        lat, lon = info.getQuerySIFTCoord(siftfile)
    stats, matchedimg, matches, combined = context.match(siftfile, matchdir, lat, lon)
    return matchedimg, matches

def draw_corr(queryimgpath, matchedimg, matches, matchoutpath=None):
    matchimgpath = os.path.join(context.dbdump, '%s.jpg' % matchedimg)
    if matchoutpath == None:
        matchoutpath = os.path.expanduser('~/client-out.jpg')
    H, inliers = corr.draw_matches(matches, queryimgpath, matchimgpath, matchoutpath)
    return H, inliers

if __name__ == '__main__':
    sift = os.path.expanduser('~/shiraz/query1/DSC_7638,37.87162,-122.27223sift.txt')
    image = os.path.expanduser('~/shiraz/query1/DSC_7638,37.87162,-122.27223.JPG')
    matchedimg, matches = match(sift)
    draw_corr(image, matchedimg, matches)
