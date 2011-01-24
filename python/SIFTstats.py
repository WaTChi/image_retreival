#! /usr/bin/env python

import os
import sys
import re
import math
#import pylab
import numpy
from optparse import OptionParser

if __name__=="__main__":
    if len(sys.argv) < 2:
        print "USAGE: {n} [options] dir".format(n=sys.argv[0])
        sys.exit()

    parser = OptionParser()
    parser.add_option("-p", "--plot", action="store_true", dest="plot",
                      default=False, help="plot histogram of # of \
                                           sift features per image.")
    parser.add_option("-f", "--file", dest="outfile", default=False,
                      help="Save plot to this file.")
    parser.add_option("-b", "--bins", dest="bins", default=10,
                      help="Set number of bins for histogram plot.")

    (options, args) = parser.parse_args()

    direc = args[0]
    if os.path.exists(direc):
        regex = re.compile(".*sift.txt$", re.IGNORECASE)
        files = os.listdir(direc)
        siftfs = [os.path.join(direc,f) for f in files if regex.match(f)]
        counts= []
        numfiles = len(siftfs)
        for fil in siftfs:
            f = open(fil, 'r')
            n = f.readline().split()[0]
            counts.append(int(n))
            f.close()

        counts = numpy.array(counts)
        mean = numpy.mean(counts)
        stddev = numpy.std(counts)
        if options.plot or options.outfile:
            (hist,bins) = numpy.histogram(counts, options.bins, new=True)
            print bins
            print hist
            step = bins[1]-bins[0]
            pylab.bar(bins[:-1], hist, width=step)
            pylab.axvline(mean, color="r")
            pylab.axvspan(mean-stddev, mean+stddev, facecolor="g", alpha=0.5)
            pylab.title("SIFT Statistics for \n{0}".format(os.path.abspath(direc)))
            pylab.xlabel("SIFT features per file.")
            if options.outfile:
                pylab.savefig(options.outfile)
            if options.plot:
                pylab.show()
                
        print "Number of SIFT files: {a}".format(a=numfiles)
        if numfiles != 0:
            print "Avg SIFT features per file: {0}".format(mean)
            print "STD Deviation of feats per file: {0}".format(stddev)
        print "Total SIFT features: {0}".format(sum(counts))
    else:
        print "Directory does not exist..."
        sys.exit(0)
