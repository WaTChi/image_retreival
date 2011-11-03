This readme outlines how to create and query against cells.

Dependencies:
  Numpy
  OpenCV
  PIL
  PyFLANN
  SciPy
  httplib

How to access git repository:
	1. add your ssh public key to the authorized keys of ericl on gorgan:
		/home/ericl/.ssh/authorized_keys
	2. clone from this url:
		$ git clone ssh://ericl@gorgan.eecs.berkeley.edu/home/ericl/query.git	

How to run python code:
NOTE: the python portion of this project is intended to be run from a linux machine
-if executable: python [executable].py [args]
-if not:  python
          import [filename]
          [filename].[methodname]([args])
          [ctrl-d to exit python interpreter when done]
    
Workflow:
  Downloading new earthmine data:

    python/cellget.py lat lon apothem outdir
      gets all earthmine bubbles in a square with length=2*apothem centered at lat, lon.
      outputs images to outdir
    
  Generating Cells:
    python/util.py/makecells lat lon length inputdir distance radius outputdir
      groups images into cells grid centered at lat,lon with specified length.
      distance: distance between neighboring cells
      radius: size of cells
      outputdir: where the cells are placed.

  Extracting SIFT Features (used for query and db):
    note: .bat scripts are intended for *Windows* machines.
    scripts/extractSIFT.bat data_directory
      make sure to check resizing option in the file. Currently, we down-sample all query images to have roughly the same resolution as our database images.
      if extracting query SIFT: make sure sizing is specified to match dataset
      if extracting dataset SIFT: make sure sizing is reasonable (depends on # of images per cell. currently we use 768x512)
    
  Querying a database:
    1. place query images in a directory
    2. extract run extractSIFT
    3. copy and modify querySystem.py:
      C.ambiguity: specifies max location ambiguity
      C.cacheEnable: should be on unless doing timing experiments
      C.cellradius: distance between neighboring cells. used for cell creation
      C.locator_function: optional pointer to location fuzzing method (currently either system.skew_location or system.load_location)
      C.maindir: specify where the project is located (defaults to /media/DATAPART2)
      C.match_callback: set to system.dump_combined_matches if outputing results for bayesian processing
      C.matchdistance: specifies max distance between two view locations capturing the same sceen
        ambiguity+matchdistance specifies the radius of the "ambiguity circle"
      C.ncells: specify max # of cells to query (should correlate w/ # of cores on machine)
      C.params.update: specfies kdtree search parameters - do not modify unless you know what you're doing
      C.resultsdir: where to put tagged outputs (~/topmatches)
      C.QUERY: specify name of query directory (assumed to be in C.maindir)
      C.topnresults: specify a list of top n results you want displayed (not needed if postprocessing w/ bayesian)
      Additional parameters can be set: see context.py/Context/init
    
    tagged output should be in C.resultsdir
    
    notes:
    -if a ground truth file is available, include it in system.check_img to get retrieval statistics
    -depending the filename/file format of your query, you may need to list your dataset as one of the special formats we handle in context.iter_queries_unfiltered
     currently we handle the following formats:
     -(default) lat,lon embeded in filename
     -lat,lon included in separate xml file (there ase typcially cell phone images generated using the Imageotag Android app)
     -no location information available (this is a debugging option intended for queries against a single cell)
    -if cells are not present, running querySystem will automatically build the neccisary index/kd-tree data structures in the local search cell if they do not exist.
    -caching is partially based on the name of your query folder/files. Be sure you don't use the same names for different query sets.
     you can clear the cach by purging the dirs/files present in Research/results
    -we're currently set up to do batch query processing. If you want query a single image, either have it be the only file in a query directory, or modify querySystem.py to call system.match instead of system.characterize
     
  Utilities:
    plotSIFT.py siftfile imagefile outfile
    plots sift features on the given image file.

    SIFTStats.py directory
    computes mean, stddev and histogram of SIFT features per file in directory.

    util.py
    contains several utility methods for manipulating cells and analysing data
    including python_to_matlab_groundTruth, which converts groundtruth files from python to matlab format.