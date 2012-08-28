import client
import system
from context import _Context, _Query
import os
from PIL import Image
import csv
import lexiconrank
import operator

class TestContext(_Context):
	def __init__(self):
		super(TestContext, self).__init__()
		self.QUERY = 'SingleImageTest'
		self.cellradius = 1
		self.dbname = 'mall_db'

	@property
	def dbdump(self):
		return self.dbname

	@property
	def dbdir(self):
		return self.dbname + '/cells-100/'

	@property
	def matchdir(self):
		p = 'test_matches'
		if not os.path.exists(p):
			os.mkdir(p)
		return p

def processSingeImage(basedir, image, cq=None):
	imagepath = basedir + '/' + image
	if cq == None:
		C = TestContext()
		C.check()
		Q = _Query()
		Q.setSensorCoord(0,0)
	else:
		C, Q = cq

	Q.jpgpath = imagepath
	Q.siftpath = os.path.splitext(imagepath)[0] + "sift.txt"
	Q.check()

	stats, matchedimg, matches, ranked, lmimg, lexdict = system.match(C, Q)
	#print ranked[0][2]
	print "Query Image: ", lexiconrank.separatePath(imagepath)[-2]
	print "Matched db image using SIFT: ", matchedimg
	print "Matched db image using lexicon: ", lmimg
	return (matchedimg, lmimg, lexdict)
	
def processMultipleImages(basedir, checkfile):
	C = TestContext()
	C.check()
	Q = _Query()
	Q.setSensorCoord(0,0)
	
	r = csv.reader(open(checkfile, 'rb'))
	imgdict = {}
	for row in r:
		if len(row) != 0:
			l = []
			for element in row:
				if len(element) == 2:
					l.append('0,0-00' + element)
				elif len(element) == 3:
					l.append('0,0-0' + element)
				else:
					l.append('0,0-' + element)
			imgdict[l[0]] = l[1:]
	
	onlyfiles = [ f for f in os.listdir(basedir) if os.path.isfile(os.path.join(basedir,f)) ]
	splitat='/.'
	
	siftscore = 0
	lmscore = 0
	combinedscore = 0
	agreescore = 0
	bestscore = 0
	total = 0
	listoffailedlm = []
	listoffailedsift  = []
	listoffailedcomb = []
	listof10 = []
	listof01 = []
	listoffailedbest = []
	for myfile in onlyfiles:
		myList = ''.join([ s if s not in splitat else ' ' for s in myfile]).split()
		if myList[-1].lower() == 'jpg' and imgdict[myList[-2]][0] != '0,0-00NA':
			siftimg, lmimg, lexdict = processSingeImage(basedir, myfile, (C,Q))
			maxsiftsore, siftdict = parseSift(myList[-2])
			#siftimg = max(siftdict.iteritems(), key=operator.itemgetter(1))[0]
			 
			total +=1
			if lmimg in imgdict[myList[-2]]:
				lmscore += 1
			else:
				listoffailedlm.append((myList[-2], lmimg))
			if siftimg in imgdict[myList[-2]]:
				siftscore += 1
			else:
				listoffailedsift.append((myList[-2], siftimg))
				
			if siftimg == lmimg:
				combimg = siftimg
			elif maxsiftsore < 10 and lmimg != None:
				combimg = lmimg
			elif lmimg != None:
				for element in siftdict.keys():
					siftdict[element] = siftdict[element] + lexdict[element]
				combimg = max(siftdict.iteritems(), key=operator.itemgetter(1))[0]
			else:
				combimg = siftimg
			
			if combimg in imgdict[myList[-2]]:
				combinedscore += 1
			else:
				listoffailedcomb.append((myList[-2], combimg, "siftimg: " + str(siftimg), siftimg in imgdict[myList[-2]], 
										"lming: " + str(lmimg), lmimg in imgdict[myList[-2]]))
				
			if siftimg in imgdict[myList[-2]] or lmimg in imgdict[myList[-2]]:
				bestscore += 1
			else:
				listoffailedbest.append(((myList[-2], lmimg, siftimg)))
				
			if lmimg in imgdict[myList[-2]] and siftimg in imgdict[myList[-2]]:
				agreescore += 1
			if lmimg in imgdict[myList[-2]] and siftimg not in imgdict[myList[-2]]:
				listof01.append((myList[-2], lmimg, siftimg))
			if siftimg in imgdict[myList[-2]] and lmimg not in imgdict[myList[-2]]:
				listof10.append((myList[-2], siftimg, lmimg))
				
	siftscore = str(siftscore)
	lmscore = str(lmscore)
	combinedscore = str(combinedscore)
	agreescore = str(agreescore)
	bestscore = str(bestscore)
	total = str(total)
	
	listoffailedlm = sorted(listoffailedlm)
	listoffailedsift = sorted(listoffailedsift)
	listoffailedcomb = sorted(listoffailedcomb)
	listof01 = sorted(listof01)
	listof10 = sorted(listof10)
	listoffailedbest = sorted(listoffailedbest)
	
	print "SiftScore: " + siftscore
	print "Lmscore: " + lmscore
	print "Combined Score: " + combinedscore
	print "Agree Score: " + agreescore
	print "Optimal Score: " + bestscore
	print "Total Num Images: " + total
	
	f = open('corytestresults.txt', 'wb')
	f.write("SiftScore: " + siftscore +'\n')
	f.write("Combined Score: " + lmscore +'\n')
	f.write("Agree Score: " + combinedscore +'\n')
	f.write("Lmscore: " + agreescore +'\n')
	f.write("Total Num Images: " + total +'\n')
	
	f.write("\nSift Failed to Find Correct Match for these Images: \n")
	for element in listoffailedsift:
		f.write(str(element) +'\n')
	f.write("\nLM Failed to Find Correct Match for these Images: \n")
	for element in listoffailedlm:
		f.write(str(element) +'\n')
	f.write("\nComb Failed to Find Correct Match for these Images: \n")
	for element in listoffailedcomb:
		f.write(str(element) +'\n')
	f.write("\nOptimal Comb Failed to Find Correct Match for these Images: \n")
	for element in listoffailedbest:
		f.write(str(element) +'\n')
	f.write("\nFound Correct Img using Sift but not Lexicon: \n")
	for element in listof10:
		f.write(str(element) +'\n')
	f.write("\nFound Correct Img using Lexicon but not Sift: \n")
	for element in listof01:
		f.write(str(element) +'\n')	
	
	f.close()
	
def parseSift(imgname):
	mydict = {}
	f = open('test_matches/' + imgname + 'sift.txt,0,0,0.res', 'rb')
	for line in f.readlines():
		mylist = line.strip().split()
		mydict[mylist[1][:-8]] = float(mylist[0])
	f.close()
	return (max(mydict.iteritems(), key=operator.itemgetter(1))[1], lexiconrank.normalizeScores(mydict))

#def parseLexicon(imgname):
#	mydict = {}
#	f = open('test_matches/' + imgname + 'lexicon.res', 'rb')
#	for line in f.readlines():
#		mylist = line.strip().split(':')
#		mydict[mylist[0].strip()] = float(mylist[1])
#	f.close()
#	total = sum(mydict.values())
#	return mydict
	
if __name__ == '__main__':
	basedir = '/home/jason/Desktop/query/project/src/tutorial/mall_query'
	checkfile = basedir + '/answers.txt'
	
	#processSingeImage(basedir, '0,0-0003.jpg')
	processMultipleImages(basedir, checkfile)
