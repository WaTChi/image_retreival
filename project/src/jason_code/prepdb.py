import sys
import os
import client
import glob
import util
from PIL import Image

def preprocesssift(dbname, imagetype):
	for path in glob.glob(dbname +'/*.'+ imagetype):
		height,width = Image.open(path).size
		
		client.preprocess_image(path, path[:-4] + '.pgm', width, height)
	for path in glob.glob(dbname +'/*.pgm'):
		client.extract_features(path)

def makecells(dbname):
	util.makecells(lat=0,lon=0,length=1,inputDir=dbname,distance=1,radius=1,outputDir=dbname+'/cells-100/')
	
def preprocesslexicon(dbname, imagetype, cannylist):
	for path in glob.glob(dbname +'/*.'+ imagetype):
		client.extract_lexicon(path, cannylist)
		
	
def fullprep(dbname, imagetype, queryname, cannylist):
#	os.system('rm -rf ' + dbname)
#	os.system('mkdir ' + dbname)
#	os.system('cp tutorial/'+dbname+'/*.jpg ' + dbname)
#	preprocesssift(dbname, imagetype)
#	preprocesssift(queryname, imagetype)
#	preprocesslexicon(dbname, imagetype, cannylist)
#	preprocesslexicon(queryname, imagetype, cannylist)
	os.system('mkdir ' + dbname + '/cells-100')
	makecells(dbname)
	
	
def lexiconprep(dbname, imagetype, queryname, cannylist, option=2):
	if option == 0 or option == 2:
		preprocesslexicon(dbname, imagetype, cannylist)
	if option == 1 or option == 2:
		preprocesslexicon(queryname, imagetype, cannylist)
	
def siftprepsingle(path):
	height,width = Image.open(path).size
	client.preprocess_image(path, path[:-4] + '.pgm', width, height)
	client.extract_features(path[:-4] + '.pgm')
	
def cleantessfiles(dir):
	os.system('rm ' + dir + '/*.tess')
	
if __name__ == '__main__':
	dbname = 'cory_pose_db'
	queryname = 'cory_pose_query'
	imagetype = 'jpg'
	fullprep(dbname, imagetype, queryname, None)

#	preprocesssift(dbname, imagetype)
#	preprocesssift(queryname, imagetype)
	
	#cleantessfiles(dbname)
	#cleantessfiles(queryname)
#	lexiconprep(dbname, imagetype, queryname, [50, 100, 150], 2)
#	client.extract_lexicon('mall_db/0,0-0003.jpg')
#	client.extract_lexicon('mall_db/0,0-0080.jpg')
	#siftprepsingle('mall_db/0,0-0081.jpg')

	
