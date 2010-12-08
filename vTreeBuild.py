#!/usr/bin/env python

import sys
import vocabularyTree
import time
import os

def buildTree(path, k=10, d=5):
#    try:
#        os.remove(os.path.join(path,"doneDatabase"+repr(k)+"_"+repr(d)+".bin"))
#    except OSError:
#        pass
    a = vocabularyTree.vocabularyTree(k, d)
    if os.path.exists(os.path.join(path,"doneSIFTTree"+repr(k)+"_"+repr(d)+".bin")):
        print "from file"
        a.fromFile(os.path.join(path,"doneSIFTTree"+repr(k)+"_"+repr(d)+".bin"))
    else:
        print "creating"
        tim = a.buildtree(path)
        print "Took {f} seconds to build the tree...{t}".format(f=tim, t=time.asctime())
        a.treeStats()
    if not os.path.exists(os.path.join(path,"doneSIFTTree"+repr(k)+"_"+repr(d)+".bin")):
        a.toFile(os.path.join(path,"doneSIFTTree"+repr(k)+"_"+repr(d)+".bin"))
    D = vocabularyTree.Database(a)
    print "inserting dir"
    D.insertdir(path)
    print"computing idf"
    D.computeIDF()
    print "normalizing"
    D.normalize()
    print "saving"
    D.toFile(os.path.join(path,"doneDatabase"+repr(k)+"_"+repr(d)+".2.bin"))
    del D
#    return D

if __name__=="__main__":
    if len(sys.argv) < 2:
        print "USAGE: {n} IMGPATH [k d]".format(n=sys.argv[0])
        sys.exit()
    try:
        dBase = buildTree(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
    except IndexError:
        dBase = buildTree(sys.argv[1])


