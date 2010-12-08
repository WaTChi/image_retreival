import vocabularyTree
import sys
import os
import re
import operator

#given a single feature, finds the best matches
#from a set of features and returns the distance
def bestMatch(feat, testFeats):
    min = float('inf')
    min2 = float('inf')
    for candidate in testFeats:
        dist = ((candidate-feat)**2).sum()
        if dist<min:
            min2=min
            min=dist
        elif dist<min2:
            min2=dist
    return float(min)/min2

# def bestMatches(feat, testFeats):
    # mins = []
    # for candidate in testFeats:
        # dist = (((candidate-feat)**2).sum())
        # mins.append(dist)
    # mins.sort()
    # return mins

#given a set of query features, returns set of best distances scores for each feature
def rankedBestMatches(queryFeats, testFeats):
    scores=[]
    for feat in queryFeats:
        scores.append(bestMatch(feat, testFeats))
    scores.sort()
    return scores

def match_scores(fi, di):
    query = vocabularyTree.readonefileSIFT(fi)
    regex = re.compile(r'.*sift.txt$', re.IGNORECASE)
    try:
        files = os.listdir(di)
    except OSError:
        raise
    files = [f for f in files if regex.match(f)]
    files.sort()
    scores={}
    for f in files:
        test = vocabularyTree.readonefileSIFT(os.path.join(di,f))
        scores[f] = rankedBestMatches(query, test)
    return scores

def add(x,y): return x+y
def firsts(x): return [x[0],x[1][0],x[1][1],x[1][2]]

def sorted_top_ratios(scores):
    results=map(firsts, scores.items())
    results.sort(key=operator.itemgetter(0))
    return results
	
def match(fi, di):
    a=match_scores(fi,di)
    r=sorted_top_ratios(a)
    print "\n".join(map(lambda x: "\t".join(map(str,x)),r))
    return r

def top_sums(scores):
    range = [1,2,5,10,15,20,25]
    sums = []
    for r in range:
        sums.append(reduce(add, scores[0:r]))
    return range, sums

if __name__=="__main__":
    if len(sys.argv) < 2:
        print "USAGE: {n} queryPath dataPath".format(n=sys.argv[0])
        sys.exit()
    match(sys.argv[1], sys.argv[2])
