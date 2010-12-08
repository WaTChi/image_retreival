#!/usr/bin/env python2.6
import sys
import os
import re
import vocabularyTree

def queryTree(tPath, qPath, outfile, numresults=10):
    if not os.path.exists(tPath):
        raise OSError
    print "Loading database...",
    d = vocabularyTree.Database(None)
    d.fromFile(tPath)
    print "Loaded."
    regex = re.compile(r'.*sift.txt$', re.IGNORECASE)
    try:
        files = os.listdir(qPath)
    except OSError:
        raise
    files = [f for f in files if regex.match(f)]
    files.sort()
    of = open(outfile, 'w')
    for f in files:
        print "Querying with {0}".format(f)
        feats = vocabularyTree.readonefileSIFT(os.path.join(qPath,f))
        results = d.queryn(feats, numresults)
        results.reverse()
        # of.write(f+","+str([r for r in results])+"\n")
        of.write(f)
        of.write(" ")
        for r in results:
            of.write(str(r[0]))
            of.write(" ")
        of.write("\n")
    of.close()

if __name__=="__main__":
    if len(sys.argv) <= 3:
        print "USAGE: {0} DBASEFILE QUERYDIR OUTFILE [numresults]".format(sys.argv[0])
        sys.exit()
        
    if len(sys.argv) == 5:
        queryTree(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        queryTree(sys.argv[1], sys.argv[2], sys.argv[3])
