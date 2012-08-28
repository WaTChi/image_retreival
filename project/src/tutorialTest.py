import client
import system
from context import _Context, _Query
import os

class TestContext(_Context):
	def __init__(self):
		super(TestContext, self).__init__()
		self.QUERY = 'SingleImageTest'
		self.cellradius = 236.6

	@property
	def dbdump(self):
		return 'testdb'

	@property
	def dbdir(self):
		return 'testdb/cells-236.6'

	@property
	def matchdir(self):
		p = 'test_matches'
		if not os.path.exists(p):
			os.mkdir(p)
		return p

image = 'tutorial/example_query_images/query1.png'
client.preprocess_image(image, image[:-4] + '.pgm', width=200, height=200)
client.extract_features(image[:-4] + '.pgm')
C = TestContext()
C.check()
Q = _Query()
Q.jpgpath = image
Q.siftpath = os.path.splitext(image)[0] + "sift.txt"
Q.setSensorCoord(37.875507, -122.264883)
Q.check()
stats, matchedimg, matches, ranked = system.match(C, Q)

print "Matched db image ", matchedimg
print "Visualization in ", C.resultsdir

