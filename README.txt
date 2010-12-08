Dependencies:
  Numpy
  SciPy
  httplib
  PIL

Important Libraries:
  earthmine.py - download data from earthmine
  vocabularyTree.py - used to construct vocab. trees and databases from earthmine image data.

Workflow:

Downloading new data:

  cellget.py lat lon radius imgx imy [depth]
    gets all earthmine bubbles within radius of lat,lon to construct a cell.
    depth information also downloaded if depth param provided.
    Download path, FOV, and images per bubble controlled by variables in cellget.py

Extracting SIFT Features:
  extractSIFT.bat data_directory
    make sure to check resizing is turned off if you don't want it. (don't want for data images)

Building Vocabulary Tree:
  vTreeBuild.py DIRECTORY k d
    Builds a vocab tree with branching factor k and depth d. Works best when k^d 1 - 2x smaller than number of features.
    Outputs two files:
      doneSIFTTreek_d.bin (not used as input for other programs)
      doneDatabasek_d.bin (used as input)

Querying a database:
  Extract SIFT features from query images (resize in extractSIFT.bat)
  THEN:
    queryDatabase.py DATABASE QUERYDIR OUTFILE [numresults]

Building a macrotree:
  Run macrotree.py
  This is a shell script. You will need to adjust parameters internally, specifically
  prefix and args, which together are the path to the four cells' SIFT files,
  as well as the path for the finished tree. (2nd to last line)
    

Querying with a macrotree:
  querymacrotree db1 db2 db3 db4 macrotree querydir resultfile

  If you want to use a different number of results, this is adjusted in querymacrotree.py on line 11


Utilities:
  plotSIFT.py siftfile imagefile outfile
    plots sift features on the given image file.

  SIFTStats.py directory
    computes mean, stddev and histogram of SIFT features per file in directory.