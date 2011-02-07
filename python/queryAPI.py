import sys
import os
import os.path
sys.path.append(os.path.expanduser("~/shiraz/Research/app/dev/python/"))
import query
import util
import corr
import client

def preprocess_image(inputfile, outputfile=None, width=768, height=512):
    """Use the convert utility to preprocess an image."""
    if outputfile == None:
        outputfile = inputfile.rsplit(".",1)[0] + ".pgm"
    os.system("convert {0} -resize {2}x{3} {1}".format(inputfile, outputfile, width, height))
	return outputfile

SIFTEXEC = os.path.join(maindir, 'Research/app/siftDemoV4/sift')

def extract_features(inputfile, outputfile=None, siftexec=SIFTEXEC):
    """Call the sift utility to extract sift features."""
    if outputfile == None:
        outputfile = inputfile.rsplit(".",1)[0] + "sift.txt"
    os.system("{0} <{1} >{2}".format(siftexec, inputfile, outputfile))
	return outputfile
