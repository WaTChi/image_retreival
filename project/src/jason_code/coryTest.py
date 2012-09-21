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
		self.dbname = None

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

def processSingeImage(dbasedir, qbasedir, image, cq=None):
	imagepath = qbasedir + '/' + image
	if cq == None:
		C = TestContext()
		C.dbname = dbasedir
		C.check()
		Q = _Query()
		Q.setSensorCoord(0,0)
	else:
		C, Q = cq

	stats, matchedimg, matches, ranked = system.match(C, Q)

	print "Query Image: ", lexiconrank.separatePath(imagepath)[-2]
	print "Matched db image using SIFT: ", matchedimg
	return matchedimg

	
if __name__ == '__main__':
	qbasedir = '/home/jason/Desktop/query/project/src/jason_code/cory_pose_query'
	dbasedir = '/home/jason/Desktop/query/project/src/jason_code/cory_pose_db'
	image = '0,0-0001.jpg'
#	checkfile = basedir + '/answers.txt'
	matchedimg = processSingeImage(dbasedir, qbasedir, image, None)
	print matchedimg
	
