#!/usr/bin/env python

#Module for building, saving, reading, searching a vocabulary tree.
#Project specific -- not for general use

#Each entry in the tree is:
# tree[x].desc an n-dimensional array that defines the vector
# tree[x].children a k-dimensional index to all the children.
#    This is 'None' when a leaf node
# tree[x].imList a path to a file containing the names of all images
#    that match this visual word


#Joe Blubaugh

import sys
import os
import re
#import kmeans
#import KLkmeans
import cPickle
import math
import time
import numpy as np
import scipy.cluster.vq as cluster
#sys.path.append("/host/earthMineData/app/lib")
sys.path.append("C:\\earthMineData\\app\\lib")
#import siftcat
import random
import struct

CHOGREGEXSTR = r'.*feate\.txt$'
SIFTREGEXSTR = r'.*sift\.txt$'
MSIFTREGEXSTR = r'.*\.msift$'
CHOGUC = "/home/blubaugh/CHoG/chog -u "
CHOGDIM = 45
KL = 0
EU = 1
MAXSIFTFEATURES = 200000
DEBUGLEVEL = 0
KMEANSITER = 25
SIFTDIM = 128

def customwhiten(data):
    #Scales each column of the data by its standard deviation.
    dprint("Whitening data...\n",1)
    (r, c) = data.shape
    dprint("Data has shape ({0}, {1})\n".format(r,c),1)
    scales = [None for i in range(c)]
    for i in range(c):
        scales[i] = data[:,i].std()
        data[:,i] /= scales[i]
        

def dprint(string, level):
    if level <= DEBUGLEVEL:
        print "DEBUG ({n}): ".format(n=level) + string+" "+time.asctime()

def uniqify(seq):
    return list(set(seq)) #Duh. Note: NOT order preserving

def featTupleCompare(tup1, tup2): #Low significance scores are 'better'
    if tup1[1] > tup2[1]:
        return -1
    elif tup1[1] == tup2[1]:
        return 0
    else:
        return 1

def readonefileSIFT(path):
    if os.path.exists(path):
        f = open(path, 'r')
        (numfeats, dim) = f.readline().split()
        numfeats = int(numfeats)
        dim = int(dim)
#        dprint("NumFeats: {n}  Dim: {d}  ".format(n=numfeats, d=dim),1)
        feats = []
        for i in range(numfeats):
            #Read first line...scale, orientation, etc
            #discard for now
            f.readline()
            count = 0
            desc = []
            while count < dim:
                l = f.readline().strip().split()
                count += len(l)
                desc.extend(l)
            try:
                feats.append(np.array(desc, dtype=np.uint8))
            except NameError:
                pass
        feats = np.array(feats, dtype=np.uint8)
        f.close()
        return feats
    else:
        raise OSError("{p} does not exist.".format(p=path))

def readonefileMSIFT(path):
    if os.path.exists(path):
        f = open(path, 'rb')
        dim = 128
        numfeats = f.read(4)
        #dprint("Raw String: {s}".format(s=numfeats),1)
        numfeats = struct.unpack("<L", numfeats)[0]
        dprint("NumFeats: {n}  Dim: {d}  ".format(n=numfeats, d=dim),1)
        feats = []
        i = 0
        while i < numfeats:
            f.read(4*4) #4 32-bit values thrown out
            d = f.read(128)
            #Conversion to 128 int array here.
            desc = [struct.unpack("B", t)[0] for t in d]
            i += 1
            try:
                feats.append(np.array(desc, dtype=np.uint8))
            except NameError:
                pass
        feats = np.array(feats, dtype=np.uint8)
        f.close()
        return feats
    else:
        raise OSError("{p} does not exist.".format(p=path))



def EUdistance(feat1, feat2):
    """Computes the Squared Euclidean Distance between two features."""
    #c = np.empty(1)
    #If they are numpy arrays:
    #if type(feat1) == type(c) and type(feat2) == type(c):
    try:
        d = feat2 - feat1
        dist = d*d
    #if not:
    except:
        dist = [(feat1[i]-feat2[i])*(feat1[i]-feat2[i]) for i in range(len(feat1))]
    return sum(dist)

def klDistance(feat1, feat2):
    """Computes KL Distance of two features."""
    c = np.empty(1)
    if type(feat1) == type(c) and type(feat2) == type(c):
        dist = feat1 * (feat1 / feat2) #All NumPy ops are element by element
    else:
        dist = [feat1[i] * math.log((feat1[i]/feat2[i]),2) for i in range(len(feat1))]
    return sum(dist)

class DocElement():
    def __init__(self, count):
        self.count = count

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "Count: "+str(self.count)

class VisualWord():
    def __init__(self, weight):
        self.weight = weight
        self.docs = {}
        self.count = 0

class Database():
    def __init__(self, tree):
        """Sets up database with tree to use for indexing."""
        self.tree = tree
        self.data = {}
        self.doccount = 0
        self.roots = {}

    def insertdir(self, path):
        """Inserts a set of feature files as individual documents."""
        ta = time.time()
        if os.path.exists(path):
            regex = re.compile(SIFTREGEXSTR, re.IGNORECASE)
            files = os.listdir(path)
            files = [f for f in files if regex.match(f)]
            readfile = readonefileSIFT
            start = 0
            stop = -8
            if len(files) == 0:
                regex = re.compile(MSIFTREGEXSTR, re.IGNORECASE)
                files = os.listdir(path)
                files = [f for f in files if regex.match(f)]
                readfile = readonefileMSIFT
                start = 5
                stop = -6
            files.sort()
            for fiel in files:
                ident = fiel[start:stop]
                feats = readfile(os.path.join(path, fiel))
                dprint("Inserting document {n}\n".format(n=ident), 1)
                self.insertdoc(feats, ident)
            tb = time.time()
            dprint("Inserting docs took {0} seconds.".format(tb-ta),1)
            ta = time.time()
            self.computeIDF()
            tb = time.time()
            dprint("Computing IDF took {0} seconds.".format(tb-ta), 1)
            ta = time.time()
            self.normalize()
            tb = time.time()
            dprint("Normalizing took {0} seconds.".format(tb-ta), 1)
        else:
            raise OSError("Path does not exist.")
        
    def insertdoc(self, feats, ident, x=None, y=None):
        """Processes a set of features as a single document, quantizes
           each feature and adds the document list to the corresponding
           visual word."""
        count = 0
        #dprint("Processing feature: {0:5}".format(count),1)
        qTime = 0
        bTime = 0
        for f in feats:
            #dprint("\b\b\b\b\b{0:5}".format(count),1)
            try:
                ta = time.time()
                vw = self.tree.quantize(f)
                tb = time.time()
            except ValueError:
                ta = time.time()
                vw = self.tree.quantize(np.transpose(f))
                tb = time.time()
            qTime += tb - ta
            ta = time.time()
            if vw in self.data:
                if ident in self.data[vw].docs:
                    self.data[vw].docs[ident].count += 1
                else:
                    self.data[vw].docs[ident] = DocElement(1)
            else:
                self.data[vw] = VisualWord(0)
                self.data[vw].docs[ident] = DocElement(1)
            tb = time.time()
            bTime += tb - ta
        self.doccount += 1
        dprint("Quantizing took {0} seconds, building structures took {1} ".format(qTime, bTime),1)

    def printdbvector(self, ident):
        """Print the vector for a particular document."""
        vector = {}
        for vw in self.data:
            if ident in self.data[vw].docs:
                vector[vw] = self.data[vw].docs[ident].count*self.data[vw].weight
        print "DB vector for {0}:".format(ident)
        print vector
        print "Normalizer: {0}".format(self.roots[ident])
        

    def computeIDF(self):
        for vw in self.data:
            self.data[vw].weight = math.log(self.doccount /
                                            float(len(self.data[vw].docs)) )

    def normalize(self):
        """Normalize Di vectors so that we can score accurately."""
        #Use L2 norm
        sums = {}
        for vw in self.data:
            for doc in self.data[vw].docs:
                try:
                    #L1
                    sums[doc] += (self.data[vw].docs[doc].count*self.data[vw].weight)
                    #L2
                    #sums[doc] += (self.data[vw].docs[doc].count * \
                    #             self.data[vw].weight) * \
                    #            (self.data[vw].docs[doc].count * \
                    #             self.data[vw].weight)
                except KeyError:
                    #L1
                    sums[doc] = (self.data[vw].docs[doc].count*self.data[vw].weight)
                    #L2:
                    #sums[doc] = (self.data[vw].docs[doc].count * \
                    #            self.data[vw].weight) * \
                    #           (self.data[vw].docs[doc].count * \
                    #            self.data[vw].weight)
        self.roots = {}
        for doc in sums:
            #L1
            self.roots[doc] = sums[doc]
            #L2
            #self.roots[doc] = math.sqrt(sums[doc])

    def queryn(self, feats, n):
        """Returns the top n matches for the doc, where doc is an iterable
           sequence of features. Results are returned in (id, significance)"""
        #Quantize each feature and collect nodes, node counts, and weights
        q = {}
        docs = {}
        sumsquared = 0
        dprint("Matching against {n} features...".format(n=len(feats)), 1)
        ta = time.time()
        #Whiten doc according to tree scales...
##        for i in range(SIFTDIM):
##            doc[i] /= self.tree.scales[i]
        for feat in feats:
            vw = self.tree.quantize(feat)
            if vw in q:
                #Increase count of path through each node
                if vw in self.data:
                    q[vw] += self.data[vw].weight
                else:
                    pass
            else:
                if vw in self.data:
                    q[vw] = self.data[vw].weight
                    docs[vw] = self.data[vw].docs
                else:
                    q[vw] = 0
                    docs[vw] = None
        #Normalize qi
        for vw in q:
            #L1
            sumsquared += q[vw]
            #L2
            #sumsquared += q[vw]*q[vw]
        #L1
        normalizer = sumsquared
        #L2
        #normalizer = math.sqrt(sumsquared)
        #Score the collected di vectors for relevance
        sig = {}
        for vw in q:
            if vw in self.data:
                for doc in self.data[vw].docs:
                    if doc in sig:
                        #L1
                        sig[doc] += (abs(q[vw] / normalizer -
                                         self.data[vw].weight *
                                         self.data[vw].docs[doc].count /
                                         self.roots[doc]) -
                                     q[vw] / normalizer -
                                     self.data[vw].weight *
                                     self.data[vw].docs[doc].count /
                                     self.roots[doc])
                        #L2
                        #sig[doc] -= 2 * (q[vw] / normalizer *
                        #                 self.data[vw].docs[doc].count *
                        #                 self.data[vw].weight /
                        #                 self.roots[doc])
                    else:
                        #L1
                        sig[doc] = 2 + (abs(q[vw] / normalizer -
                                       self.data[vw].weight *
                                       self.data[vw].docs[doc].count /
                                       self.roots[doc]) -
                                   q[vw] / normalizer -
                                   self.data[vw].weight *
                                   self.data[vw].docs[doc].count /
                                   self.roots[doc])
                        #L2
                        #sig[doc] = 2 - \
                        #           2 * (q[vw] / normalizer *
                        #                self.data[vw].docs[doc].count *
                        #                self.data[vw].weight /
                        #                self.roots[doc])
        

        #Turn sig into a list
        sig = sig.items()
        sig.sort(featTupleCompare)
        #Return the top N
        tb = time.time()
        dprint("All matching took {t} seconds.\n".format(t=(tb-ta)), 1)
        return sig[(-1*n):]
    

    def toFile(self, path):
        f = open(path, 'wb')
        self.tree.toFile(path+".tree")
#        cPickle.dump(self.tree,f)
        cPickle.dump(self.data,f)
        cPickle.dump(self.roots,f)
        cPickle.dump(self.doccount,f)
        f.close()

    def fromFile(self, path):
        f = open(path, 'r')
        self.tree = vocabularyTree()
        self.tree.fromFile(path+".tree")
#        self.tree = cPickle.load(f)
        self.data = cPickle.load(f)
        self.roots = cPickle.load(f)
        self.doccount = cPickle.load(f)
        f.close()

    def dbStats(self):
        self.tree.treeStats()
        print "Visual words in database: {0}".format(len(self.data))

class vTreeNode():
    def __init__(self, desc=None, children = None):
        #If a node has an imList, it is a leaf node
        self.desc = desc
        self.children = children #should be a list of vTreeNodes
        self._entropy = 0 #Should only me modified when entropy is
                          #computed for the tree
        self.index = -1

            

class vocabularyTree():
    def __init__(self, k=0, D=0):
        self.tree = []
        #A kD tree can be represented in a flat array using indices:
        #    childI = parentI x k + {1, 2, ..., k}
        if type(k) is int:
            self.k = k
        if type(D) is int:
            self.D = D
        self.dMeasure = EU
        self.rawData = []
        self.scales = None
        self.chogregex = re.compile(CHOGREGEXSTR)
        self.siftregex = re.compile(SIFTREGEXSTR)
        self.msiftregex = re.compile(MSIFTREGEXSTR)

    def readRawSIFTData(self, path):
        """Reads in a SIFT feature set from XYZ"""
        self.rawData = []
        self.dMeasure = EU
        try:
            files = os.listdir(path)
            os.chdir(path)
        except OSError:
            raise

        feats = [f for f in files if self.siftregex.match(f)]
        readfile = readonefileSIFT
        
        if len(feats) == 0: #Get MSIFT feats
            feats = [f for f in files if self.msiftregex.match(f)]
            readfile = readonefileMSIFT 
        del files
        for x in feats:
            data = readfile(x)
            self.rawData.extend(data)
#        self.rawData = np.array(self.rawData)
        dprint("Done reading in features...\n",1)

    def addRawSIFTData(self, path):
        #Don't clear raw data before starting.
        rD = []
        try:
            files = os.listdir(path)
            os.chdir(path)
        except OSError:
            raise

        feats = [f for f in files if self.siftregex.match(f)]
        readfile = readonefileSIFT
        
        if len(feats) == 0: #Get MSIFT feats
            feats = [f for f in files if self.msiftregex.match(f)]
            readfile = readonefileMSIFT 
        del files
        for x in feats:
            data = readfile(x)
            rD.extend(data)
        rD = np.array(rD)
        self.rawData = np.concatenate((self.rawData, rD), 0)
        dprint("Done reading in features...\n",1)

    def buildtree(self, datapath, k="k", D="D"):
        """Performs k-means clustering D times on self.rawData
           to produce a vocabulary tree."""
        if k is "k":
            k = self.k
        else:
            self.k = k

        if D is "D":
            D = self.D
        else:
            self.D = D
        t1 = time.time()
        if os.path.exists(datapath):
            if len(self.rawData) == 0:
                self.readRawSIFTData(datapath)
            if len(self.rawData) > MAXSIFTFEATURES:
                indices = random.sample(range(len(self.rawData)), MAXSIFTFEATURES)
                self.tData = [self.rawData[i] for i in indices]
                self.rawData = np.array(self.tData)
            else:
                self.rawData = np.array(self.rawData)
            #Whiten data
            self.scales = customwhiten(self.rawData)
            self.tree = self._buildtree(self.rawData, None, self.k, self.D)
            self.rawData = None
            t2 = time.time()
            self.treeStats()
            return t2 - t1
        return -1
        
    def _buildtree(self, data, desc, k, D):
        """Performs k-means clustering D times on self.rawData
           to produce a vocabulary tree."""
        
        if (D == 0) or (len(data) < k): #Done recursing
            return vTreeNode(desc = desc, children = None)
        else:
            dprint("Running kmeans with {X} points...".format(X=len(data)),1)
            try:
                (centroids, labels) = cluster.kmeans2(data, k, minit="points",iter=KMEANSITER)
                #centroids = cluster.kmeans(data, k, iter=KMEANSITER)
                #labels = cluster.vq(data, centroids)
            except UserWarning: #Repeat once if we get the 'empty cluster' warning
                (centroids, labels) = cluster.kmeans2(data, k, minit="points",iter=KMEANSITER)
                #centroids = cluster.kmeans(data, k, iter=KMEANSITER)
                #labels = cluster.vq(data, centroids)
            childList = []
            dStore = [np.empty(0) for i in range(k)]
            for i in range(len(data)):
                dStore[labels[i]] = np.append(dStore[labels[i]], data[i], axis=0)

            #reshape and recurse
            for i in range(k):
                dStore[i] = np.reshape(dStore[i], (len(dStore[i])/SIFTDIM,SIFTDIM))
                childList.append(self._buildtree(dStore[i], centroids[i], k, D-1))
            return vTreeNode(desc = desc, children = childList)        
        
    def toFile(self, path):
        """Writes vocabulary tree to a binary file descriptor."""
        f = open(path, "wb")
        cPickle.dump(self.k,f)
        cPickle.dump(self.D,f)
        cPickle.dump(self.dMeasure, f)
        cPickle.dump(self.numZeros, f)
        cPickle.dump(self.numLeaf, f)
        cPickle.dump(self.numNodes, f)
        cPickle.dump(self.scales, f)
        try:
            cPickle.dump(self.tree, f)
        except NameError:
            pass
        f.close()

    def fromFile(self, path):
        """Reads vocabulary tree from a binary file, replacing
           current values."""
        try:
            f = open(path, 'r')
        except IOError:
            print "Could not open file {F} for reading...".format(F = path)
            raise
        else:
            self.k = cPickle.load(f)
            self.D = cPickle.load(f)
            self.dMeasure = cPickle.load(f)
            self.numZeros = cPickle.load(f)
            self.numLeaf = cPickle.load(f)
            self.numNodes = cPickle.load(f)
            self.scales = cPickle.load(f)
            self.tree = cPickle.load(f)
            #Parse out object list
            f.close()


    def quantize(self, feature):
        """Returns an integer that quantizes the features path down the tree.
           Feature must have the same dimensionality as the tree."""
        nearest = self.tree
        depthLeft = self.D
        shiftamt = int(math.ceil(math.log(self.k, 2)))
        distance = EUdistance
##        try:
##            if self.dMeasure == KL:
##                distance = kldistance
##            else:
##                distance = EUdistance
##        except NameError:
##            distance = EUdistance
        quant = 0
        while nearest.children is not None: #want leaf node
            dist = sys.float_info[0] #max system float
            for (ind, child) in enumerate(nearest.children):
                tDist = distance(feature, child.desc)
                if tDist < dist:
                    dist = tDist
                    nNearest = child
                    index = ind
            nearest = nNearest
            depthLeft -= 1
            quant = (quant << shiftamt) + index

        if depthLeft > 0:
            quant <<= shiftamt*depthLeft

        return quant

    def getNode(self, ident):
        """Return the node with the path quantized by ident"""
        depthLeft = self.D
        shiftamt = int(math.ceil(math.log(self.k, 2)))
        mask = 1
        for i in range(0, shiftamt-1):
            mask = (mask << 1) | 1

        node = self.tree
        while depthLeft > 0:
            add = (ident >> (shiftamt * depthLeft-1)) & mask
            if node.children is not None:
                node = node.children[add]
            else:
                return node
            depthLeft -= 1

        return node
            

    def treeStats(self):
        """Reports some things about the tree, including:
           number of nodes
           number of leaf nodes
           number of nodes with '0' descriptors"""

        (self.numNodes, self.numLeaf, self.numZeros) = \
                        self._treeStats(self.tree)
        #Print results
        print "Statistics for tree:"
        print "Num Nodes:      {n}".format(n = self.numNodes)
        print "Num Leaf Nodes: {n}".format(n = self.numLeaf)
        print "Num Zeros:      {n}".format(n = self.numZeros)

    def _treeStats(self, node):
        numNodes = 1
        numLeaf = 0
        numZeros = 0
        if node.children is None:
            numLeaf = 1

        if node.desc is None: #top Node
            pass
        elif sum(node.desc) < 0.0000001:
            numZeros = 1
        else:
            pass
            
        try:
            for child in node.children:
                (tNode, tLeaf, tZero) = self._treeStats(child)
                numNodes += tNode
                numLeaf += tLeaf
                numZeros += tZero
        except TypeError:
            pass

        return (numNodes, numLeaf, numZeros)
    


if __name__=="__main__":
    #Run Test code
    pass
