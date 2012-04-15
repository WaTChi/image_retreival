This readme outlines how to create and query against cells.

Dependencies:
  Numpy
  OpenCV
  PIL
  pyFLANN
  SciPy
  httplib

On a Ubuntu machine, most of the dependencies can be installed using this command:

	$ sudo apt-get install python-numpy python-opencv python-imaging python-scipy python-httplib2

You will also have to install pyFLANN, which can be found here:

	http://www.cs.ubc.ca/~mariusm/index.php/FLANN/FLANN

How to run python code:
NOTE: the python portion of this project is intended to be run from a linux machine
-if executable: python [executable].py [args]
-if not:  python
          import [filename]
          [filename].[methodname]([args])
          [ctrl-d to exit python interpreter when done]
    
To get started, see the mini-tutorial in:

	project/src/tutorial

for information on how to setup a basic image query pipeline from scratch.
