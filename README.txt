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
  Downloading new earthmine data:

    python/cellget.py lat lon radius outdir
      gets all earthmine bubbles in a square with length=2*radius centered at lat, lon.
      outputs images to outdir
    
  Generating Cells:
    python/util.py/makecells lat lon length inputdir distance radius outputdir
      groups images into cells grid centered at lat,lon with specified length.
      distance: distance between neighboring cells
      radius: size of cells
      outputdir: where the cells are placed.

  Extracting SIFT Features:
    scripts/extractSIFT.bat data_directory
      make sure to check resizing option in the file.
    
  Querying a database:
    1. place query images in a directory
    2. extract run extractSIFT
    3. copy and modify querySystem.py:
      C.QUERY: specify location of query directory
      C.params.update: specfies kdtree search parameters - do not change unless you know what you're doing
      C.ambiguity: specifies max location ambiguity
      C.matchdistance: specifies max distance between two view locations capturing the same sceen
        ambiguity+matchdistance specifies the radius of the "ambiguity circle"
      C.topnresults: specify a list of top n results you want displayed (not needed if postprocessing w/ bayesian)
      C.ncells: specify max # of cells to query (should correlate w/ # of cores on machine)
      C.cacheEnable: should be on unless doing timing experiments
      C.location_function: optional pointer to location fuzzing method (currently either system.skew_location or system.load_location)
      C.match_callback: set to system.dump_combined_matches if outputing results for bayesian processing
      Additional parameters can be set: see context.py/Context/init
    
    tagged output should be in ~/topmatches

  Utilities:
    plotSIFT.py siftfile imagefile outfile
    plots sift features on the given image file.

      SIFTStats.py directory
      computes mean, stddev and histogram of SIFT features per file in directory.
