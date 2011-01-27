Dependencies:
  Numpy
  OpenCV
  PIL
  PyFLANN
  SciPy
  httplib

Important Libraries:
  earthmine.py - download data from earthmine

Workflow:

Downloading new data:

  cellget.py lat lon radius imgx imy [depth]
    gets all earthmine bubbles within radius of lat,lon to construct a cell.
    depth information also downloaded if depth param provided.
    Download path, FOV, and images per bubble controlled by variables in cellget.py

Extracting SIFT Features:
  extractSIFT.bat data_directory
    make sure to check resizing is turned off if you don't want it. (don't want for data images)

Querying a database:
  Extract SIFT features from query images (resize in extractSIFT.bat)
  THEN:
    querySystem.py [args]

Utilities:
  plotSIFT.py siftfile imagefile outfile
    plots sift features on the given image file.

  SIFTStats.py directory
    computes mean, stddev and histogram of SIFT features per file in directory.
