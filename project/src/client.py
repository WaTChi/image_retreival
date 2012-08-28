# There are some nice utility functions for preprocessing an image here.
#
# However, this script is deprecated.
#
# You should assemble a single-image pipeline
# using context.py and pieces from system.py
#
# I think singleImage.py is an example of a single-image context.

import os
import corr
import info
import time
from config import *
import query
#import queryContext as context

#matchdir = '/tmp/client/%s' % query.searchtype(context.params)
#if not os.path.exists(matchdir):
#    os.makedirs(matchdir)
SIFTEXEC = os.path.join('/home/jason/Desktop/query/project/src/siftDemoV4/sift')
#context.vars_init()

def preprocess_image(inputfile, outputfile=None, width=768, height=512, fill=False):
    """
    Use the convert utility to preprocess an image.
    The fill argument specifies whether or not the image should take up the full dimension
    specified. If so, cropping may occur. Default behavior is for the image to be resized
    just small enough that it fits in a width x height box.
    """
    timer = time.time()
    if outputfile == None:
        outputfile = inputfile.rsplit(".",1)[0] + ".pgm"
    if fill:
        fillparams = "^ -gravity center -extent {0}x{1}".format(width, height)
    else:
        fillparams = ""
    os.system("convert {0} -resize {2}x{3}{4} {1}".format(inputfile, outputfile, width, height, fillparams))
    INFO("--> Image conversion: " + str(time.time()-timer) + "s")
    return outputfile

def extract_lexicon(inputfile, thresholds=[130], outputfile=None):
    timer = time.time()
    imagename = inputfile.rsplit(".",1)[0]
    realname = imagename.split('/')[-1]
#    imagepath = imagename.split('/')[:-1]
    if outputfile == None:
        outputfile = imagename + ".tess"
    print outputfile
    if os.path.lexists(outputfile) == True:
        return outputfile
    maxsize = -1
    for threshold in thresholds:
        temp = "./read_text " + inputfile + ' 1 ' + str(threshold) + ' 0'
        print temp
        os.system(temp)
        if os.path.getsize('detected_text.txt') > maxsize:
            maxsize = os.path.getsize('detected_text.txt')
            os.system("cp detected_text.txt final_text.txt")
    os.system("mkdir " + imagename)
    os.system("mv correctedpatch* patch* result" + realname + ".jpg " + imagename)
    os.system("mv final_text.txt " + imagename + ".tess")
    INFO("--> Lexicon extraction: " + str(time.time()-timer) + "s")
    return outputfile

def extract_features(inputfile, outputfile=None, siftexec=SIFTEXEC):
    """Call the sift utility to extract sift features."""
    timer = time.time()
    if outputfile == None:
        outputfile = inputfile.rsplit(".",1)[0] + "sift.txt"
    if os.path.lexists(outputfile) == True:
        return outputfile
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
