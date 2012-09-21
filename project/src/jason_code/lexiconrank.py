from os import listdir
from os.path import isfile, join
import re
import os
import difflib
import operator

def addImagetoList(listofimages, cellpath):
    onlyfiles = [ f for f in listdir(cellpath) if isfile(join(cellpath,f)) ]
    for file in onlyfiles:
        myList = separatePath(file)
        if myList[-1].lower() == 'jpg' and myList[-2] not in listofimages:
            listofimages.append(myList[-2])
    return listofimages

def listallimages():
    return addImagetoList([], 'corydb')

def separatePath(path):
    splitat='/.'
    myList = ''.join([ s if s not in splitat else ' ' for s in path]).split()
    return myList

def readTessFile(imgpath, conf=None):
    imgpath = imgpath + '.tess'
    f = open(imgpath, 'rb')
    mystring = ''
    for line in f.readlines():
        mystring += re.sub('\s', '', line)
    f.close()
    if len(mystring) == 0:
        os.system('echo ' + '"' + imgpath + '" >> listofnooutput')
    return mystring

def readTessFile2_old(imgpath, conf=None):
    imgpath = imgpath + '.tess'
    f = open(imgpath, 'rb')
    mystring = []
    for line in f.readlines():
        mystring.append(re.sub('\s', '', line))
    f.close()
    if len(mystring) == 0:
        os.system('echo ' + '"' + imgpath + '" >> listofnooutput')
    return mystring

def readTessFile2(imgpath, conf=None):
    imgpath = imgpath + '.tess'
    f = open(imgpath, 'rb')
    mystring = []
    tempstring = ""
    for line in f.readlines():
        temp = re.sub('\s', '', line)
        pat = re.compile("##\d\d")
        testconf = re.search(pat, temp)
        
        if testconf:
            if conf != None and int(testconf.string.strip('#')) >= conf and len(tempstring) != 0:
                mystring.append(tempstring)
            elif conf == None and len(tempstring) != 0:
                mystring.append(tempstring)
            tempstring = ""
        else:
            tempstring += temp
    f.close()
    return mystring


def computeRanks(dbbase, listofimages, qimgPath):
    scoredict = {}
    for i, image in enumerate(listofimages):
        imagepath = dbbase+ '/' +image
        dstring = readTessFile(imagepath, 50)
        
        if i == 0:
            qstring = readTessFile(qimgPath[:-4], 50)
        
        if qstring == '':
            return None
        
        scoredict[image] = scoreMethod(dstring, qstring)
        
    myList = separatePath(qimgPath)
    f = open('test_matches/' + myList[-2] + 'lexicon.res', 'wb')
    sortlist = sorted(scoredict.items())
    for key, value in sortlist:
        f.write(str(key) + " : " + str(value)+ '\n')
    f.close()
    return scoredict

def compareTwoImages(img1, img2):
    dstring = readTessFile2(img2, 0);
    qstring = readTessFile2(img1, 0)
        
    s = scoreMethod(dstring, qstring)
    
        
    print "query string:"
    print qstring
    
    print "database string:"
    print dstring
    
    print "Similarity Ratio: " + str(s)
    
def scoreMethod(x,y):
#    for i,e in enumerate(x):
#        x[i] = re.sub(r'\W+', '', e)
#    for i,e in enumerate(y):
#        y[i] = re.sub(r'\W+', '', e)
        
    x = re.sub(r'\W+', '', x)
    y = re.sub(r'\W+', '', y)
    
    return similarityScore_trigram(x,y)
    
def similarityScore_difflib(a, b):
    if len(a) > len(b):
        return similarityScore_difflib(b, a)
    if len(a) == 0:
        return 0.0
    x = len(a)
    y = len(b)
    s = difflib.SequenceMatcher()
    maxscore = 0.0
    for i in range(0, y - x):
        s.set_seqs(a, b[i:i+x])
        if s.ratio() > maxscore:
            maxscore = s.ratio()
    return maxscore

def similarityScore_difflib_simple(x, y):
    s = difflib.SequenceMatcher(a=x,b=y)
    return s.ratio()

def similarityScore_combined(x, y):
    return similarityScore_editdist(x,y) + similarityScore_trigram(x,y)

def similarityScore_trigram(a,b):
    def returnTrigrams(s):
        l = []
        s = "  " + s + "  "
        for i in range(0, len(s) - 3):
            l.append(s[i:i+3])
        return l
    def returnTrigrams2(s):
        l = []
        for mystring in s:
            mystring = "  " + mystring + "  "
            for i in range(0, len(mystring) - 3):
                l.append(mystring[i:i+3])
        return l  
    x = returnTrigrams(a)
    y = returnTrigrams(b)
    score = 0
    for element in x:
        if element in y:
            #print element
            y.remove(element)
            score += 1
    return score
#    print set(x)
#    print set(y)
#    return len(set(x) & set(y))

def similarityScore_editdist(s1, s2):
    l1 = len(s1)
    l2 = len(s2)

    matrix = [range(l1 + 1)] * (l2 + 1)
    for zz in range(l2 + 1):
        matrix[zz] = range(zz,zz + l1 + 1)
    for zz in range(0,l2):
        for sz in range(0,l1):
            if s1[sz] == s2[zz]:
                matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz])
            else:
                matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz] + 1)
                
    return float(1.0/(matrix[l2][l1]+0.0001))
    
def normalizeScores(scoredict):
    totalscore = float(sum(scoredict.values()))
    if totalscore == 0:
        return scoredict
    
    for key in scoredict.keys():
        scoredict[key] = scoredict[key]/totalscore
    return scoredict

    
def returnTopMatch_random(dbbase, listofimages, qimgPath):
    myDict = computeRanks(dbbase, listofimages, qimgPath)
    if myDict == None:
        return (None, None)
    topimage = max(myDict.iteritems(), key=operator.itemgetter(1))[0]
    topscore = max(myDict.iteritems(), key=operator.itemgetter(1))[1]
    if topscore <= 1.0:
        return (None, None)
    return (normalizeScores(myDict), topimage)

        
if __name__ == '__main__':
    img1 = '/home/jason/Desktop/query/project/src/tutorial/cory_2nd_nexus/0,0-0078'
    img2 = '/home/jason/Desktop/query/project/src/corydb2/0,0-0064'
    dbbase = 'corydb'
    #print returnTopMatch_random(dbbase, listallimages(), img1 + '.jpg')
    compareTwoImages(img1, img2)