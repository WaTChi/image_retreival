#!/usr/bin/env python

import sys
sys.path.append("C:\\earthMineData\\app\\lib")
import vocabularyTree
sys.path.pop()
import os
import re
import time

NUMRESULTS = 15

def queryMacTree(d1, d2, d3, d4, mac, qdir, rfile):
    ts = time.time()
    path1 = d1[0:d1.rindex("doneData")]
    path2 = d2[0:d2.rindex("doneData")]
    path3 = d3[0:d3.rindex("doneData")]
    path4 = d4[0:d4.rindex("doneData")]
    D1 = vocabularyTree.Database(None)
    D2 = vocabularyTree.Database(None)
    D3 = vocabularyTree.Database(None)
    D4 = vocabularyTree.Database(None)
    ta = time.time()
    print "Reading db {0} from file...".format(d1),
    D1.fromFile(d1)
    tb = time.time()
    print "done. Took {0} seconds.".format(tb-ta)
    print "Reading db {0} from file...".format(d2),
    D2.fromFile(d2)
    ta = time.time()
    print "done. Took {0} seconds.".format(ta-tb)
    print "Reading db {0} from file...".format(d3),
    D3.fromFile(d3)
    tb = time.time()
    print "done. Took {0} seconds.".format(tb-ta)
    print "Reading db {0} from file...".format(d4),
    D4.fromFile(d4)
    ta = time.time()
    print "done. Took {0} seconds.".format(ta-tb)
    Mac = vocabularyTree.vocabularyTree()
    print "Reading macrotree {0} from file...".format(mac),
    Mac.fromFile(mac)
    tb = time.time()
    print "done. Took {0} seconds.".format(tb-ta)

    #For each sift or msift file in the query directory, query all four trees,
    #compile results, build a temp. database, get new top results.
    sregex = re.compile(vocabularyTree.SIFTREGEXSTR, re.IGNORECASE)
    msregex = re.compile(vocabularyTree.MSIFTREGEXSTR, re.IGNORECASE)
    files = os.listdir(qdir)
    fles = [f for f in files if sregex.match(f)]
    readfile = vocabularyTree.readonefileSIFT
    suffix = "sift.txt"
    start = 5
    stop = -8
    if len(fles) == 0:
        fles = [f for f in files if msregex.match(f)]
        readfile = vocabularyTree.readonefileMSIFT
        suffix = ".msift"
        stop = -6
    of = open(rfile, "w")
    for f in fles:
        print "Querying with {0}".format(f)
        ta = time.time()
        feats = readfile(os.path.join(qdir, f))
        r1 = D1.queryn(feats, NUMRESULTS)
        r1.reverse()
        #Check for cell 1 matches here...
        iden = int(f[start:stop])
        of.write(f+","+d1)
        for r in r1:
            of.write(",")
            of.write(str(r[0]))
        of.write("\n")
        r2 = D2.queryn(feats, NUMRESULTS)
        r2.reverse()
        of.write(f+","+d2+","+str([r[0] for r in r2])+"\n")
        r3 = D3.queryn(feats, NUMRESULTS)
        r3.reverse()
        of.write(f+","+d3+","+str([r[0] for r in r3])+"\n")
        r4 = D4.queryn(feats, NUMRESULTS)
        r4.reverse()
        of.write(f+","+d4+","+str([r[0] for r in r4])+"\n")
        tb = time.time()
        print "Initial queries took {0} seconds.".format(tb-ta)
        MD = vocabularyTree.Database(Mac)
        #Insert results into MD...
        #Assume we're using the same SIFT file format as the queries
        ta = time.time()
        for r in r1:
            doc = readfile(os.path.join(path1, "image{0:0>4}".format(r[0])+suffix))
            ident = "1_"+str(r[0])
            MD.insertdoc(doc, ident)
        for r in r2:
            doc = readfile(os.path.join(path2, "image{0:0>4}".format(r[0])+suffix))
            ident = "2_"+str(r[0])
            MD.insertdoc(doc, ident)
        for r in r3:
            doc = readfile(os.path.join(path3, "image{0:0>4}".format(r[0])+suffix))
            ident = "3_"+str(r[0])
            MD.insertdoc(doc, ident)
        for r in r4:
            doc = readfile(os.path.join(path4, "image{0:0>4}".format(r[0])+suffix))
            ident = "4_"+str(r[0])
            MD.insertdoc(doc, ident)
        MD.computeIDF()
        MD.normalize()
        results = MD.queryn(feats, NUMRESULTS)
        results.reverse()
        #Check again for matches...
        of.write(f+",macro")
        for r in results:
            of.write(str(r[0]))
        of.write("\n")
        tb = time.time()
        print "Insertion & queries took {0} seconds.".format(tb - ta)
    of.close()
    

if __name__=="__main__":
    if len(sys.argv) != 8:
        print "USAGE: {0} db1 db2 db3 db4 macrotree querydir resultsfile".format(sys.argv[0])
        sys.exit()
    queryMacTree(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
                 sys.argv[6], sys.argv[7])
    
