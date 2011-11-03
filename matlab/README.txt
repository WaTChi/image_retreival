This readme outlines how to use querySystem.py results to run Bayes post-processing in matlab on Windows.

Variables:
	<dir>	: String describing the root directory for most matlab files, <dir> = 'Z:\Research\app\code\matlab'
	<set>	: String describing the query set used in the run ( e.g. 'query1' )
	<dec>	: String describing the post-processing decision used ( e.g. 'bayes' )
	<dist>	: String describing the simulated ambiguity distribution ( e.g. 'exact' )
	<cld>	: The cell distance used for determining whether or not to search a cell ( e.g. 336.6 meters )
	<scld>	: String describing the cell distance <cld> equal to the nearest integer in meters ( e.g. '337' )
	
	
	
	<cdiv>	: Vector dividing the query set into divisions based on the number of candidates ( more on this later )
				e.g. <cdiv> = [10 100] to split queries into 10- , 10-100 , 100+ candidates, <cdiv> = [] for no divisions
	<scdiv>	: String describing the candidate division vector <cdiv>, e.g.
				if <cdiv> = [10 100], <scdiv> = '10,100
	

Prerequisite actions:
	Map \\gorgan\data to your Z drive
		1) Right-click on Computer in Windows Explorer; select 'Map Network Drive...'
		2) Choose drive Z: and enter \\gorgan\data as the folder
		3) Connect using your local account on gorgan ( e.g. gorgan\username )
	Copy over .res results files for your set from gorgan to your local machine's C drive
		1) See python README if you do not know where the results files have been stored on gorgan
		2) Make sure you are copying the cell combination specific results files, not the results files for each cell
		3) Copy results files and move them to C:\matlab_local\results\<set>
	Copy over .txt ground truth file for your set
		1) See python README for generating .txt ground truth files from .py ground truth files
		2) Copy your set ground truth file to C:\matlab_local\ground-truth
		3) Rename your ground truth file to gt_<set>.txt

Directory:
	Z:\Research\app\code\matlab ...
		\bayes			 	: bayes framework and functions
		\post-processing 	: post-processing framework and functions
			\trainQueryClassifier.m	 ( more on this later )
			\post_process.m	  		 ( more on this later )
		\util 			 	: utility functions
	
Main functions:
	trainQueryClassifier.m	: Main function for training the Bayes classifier
		(inp) method	 	: Structure containing the run parameters  ( more on this later )
		(inp) reset		 	: Reset boolean ( 1 = erase previous classifier , 0 = add to previous classifier - default )
		(out) classifier 	: Structure containing the classifier data from training  ( more on this later )
	post_process.m			: Main function for post-processing with Bayes classifier
		(inp) method	 	: Structure containing the run parameters  ( more on this later )
		(inp) reset		 	: Reset boolean ( 1 = erase previous results , 0 = add to previous results - default )
		(out) results	 	: Structure containing the results data from post-processing  ( more on this later )
	
(INPUT)	
method: Structure with the following contents...
	.set					: String determining which set to use for the run ( identical to the <set> variable )
	.cell_dist				: Maximum distance between location and cell center to search a cell
								This should be ambiguity radius + cell radius + 25  ( e.g. 75 + 236.6 + 25 = 336.6 )
								This is identical to the variable <cd> listed at top
								This is used to determine the string <scd>  ( e.g. <cd> = 336.6 means <scd> = '337' )
	.decision				: String in the format '<dec>-<prm>'
								<dec> is a string indicating the decision used  ( same as <dec> variable listed at top )
									'vote' ranks candidates based on vote only
									'bayes' uses the Bayes classifier to rank the candidates
								<prm> is a list of parameters used for the 'bayes' decision
									'd' indicates only distances are used for classification
									'v' indicates only votes are used for classification
									'dv' indicates to use both distances and votes for classifiction
								In most cases use method.decision = 'bayes-dv' for the full post-processing
								Note that this string is not used or necessary for trainQueryClassifier.m
	.canddiv				: Vector to train separately for for queries with different numbers of candidates
								e.g. method.canddiv = [10 100] trains separately for the following three types of queries
									those with 10 or fewer candidates
									those with 100 or fewer and more than 10 candidates
									those with more than 100 candidates
							  In most cases use method.canddiv = [] for training and post-processing ( no separation )
	.distribution			: String in the format '<dist>-<prm>' used to simulate location ambiguity
								<dist> is a string indicating the simulated distribution  ( same as <dist> listed at top )
									'exact' indicates no ambiguity simulation and has no parameter value ( e.g. 'exact-' )
									'unif' indicates a uniform distribution with parameter indicating max ambiguity
									'expo' indicates an exponential distribution with an exponential parameter
								<prm> indicates the parameter value for the simulated distribution given
								Common inputs
									'exact-'  : most frequently used ( use when you don't want to simulate ambiguity )
									'unif-75' : use to simulate uniform ambiguity up to 75 meters away
									'expo-50' : use to simulate ambiguity with an exponential distribution ( 50m param )


To run post_process.m:
	1) Open Matlab.
	2) Navigate the current directory to Z:\Research\app\code\matlab\post-processing
	3) Set the method variable and reset variable, e.g.
			method.set = 'query1';
			method.cell_dist = 336.6;
			method.canddiv = [];
			method.decision = 'bayes-dv';
			method.distribution = 'exact-';
			reset = false;
	4) Run the function, e.g.
			results = post_process(method,reset);
	The output 'results' is automatically saved to disk in <dir>\post-processing\<dec>\<set>_<dist><scd>results.mat
		e.g. above results saved in Z:\Research\app\code\matlab\post-processing\bayes\query1_exact337results.mat
		More on the results structure listed later

To run trainQueryClassifier.m:
	1) Open Matlab.
	2) Navigate the current directory to Z:\Research\app\code\matlab\post-processing
	3) Set the method variable and reset variable, e.g.
			method.set = 'query2';
			method.cell_dist = 336.6;
			method.canddiv = [];
			method.distribution = 'exact-';
			reset = true;
	4) Run the function, e.g.
			classifier = trainQueryClassifier(method,reset);
	The output 'classifier' is automatically saved to disk in 	
		<dir>\post-processing\bayes\classifier\<set><scd>_<dist><scd>results.mat
		e.g. above results saved in Z:\Research\app\code\matlab\post-processing\bayes\query1_exact337results.mat
		More on the results structure listed later

Dependencies:
  Numpy
  OpenCV
  PIL
  PyFLANN
  SciPy
  httplib

How to run matlab code:
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
    -we're currently set up to do batch query processing. If you want query a single image, set C.selection = ['substringtomatch']
     
  Utilities:
    plotSIFT.py siftfile imagefile outfile
    plots sift features on the given image file.

    SIFTStats.py directory
    computes mean, stddev and histogram of SIFT features per file in directory.

    util.py
    contains several utility methods for manipulating cells and analysing data
    including python_to_matlab_groundTruth, which converts groundtruth files from python to matlab format.
