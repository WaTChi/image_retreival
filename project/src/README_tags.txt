== Overview ==

We project tags from 3d space onto database images. Then, we move them onto query images using the computed homography.

The entry point to running a query is [querySystem.py] or [querySystemCopy.py] (there are two just for the sake of convenience). This imports the Context object, sets the parameters, and runs the query using [system.py].

[system.py] can be thought of as collection of functions that operate on individual queries (and the characterize() function calls these in the right order).

characterize() calls match(C, Q) for each individual query image.

match()	runs a query [query.py] and gets back a list of feature matches.

These feature matches are put through a series of strong filtering stages to discard bad feature matches, and are reranked.

Then the top image is picked. This would be the "topN" result.

Following that, homographies are computed between the query and database matches. This process continues down the topN list until a "good" homography match is determined. Then tags are transferred onto this image and returned to the client.

Technical details below:

== Disk Storage ==

The main data storage format used is the numpy record array, which packs structured data on disk in an efficient format. These arrays have types called dtypes, and can be queried by names inside their dtype as well as indexed names. For example, say you had the dtype of an array of 3d points in space:

    map3d_dtype = {
      'names': ['lat', 'lon', 'alt'],
      'formats': ['float64', 'float64', 'float64'],
    }

This means that each row of the record array will have three 64 bit floats. If you wanted to access the latitude third item in the array you would say:

	array[3]['lat'] or array['lat'][3]

(see numpy record array documentation for advanced operations).

If general, you want to avoid operating on record arrays using python constructions like for loops. This is very slow. Instead, try to find a numpy operation that captures what you are trying to do (select a subset, reshape the array, etc).

== What uses numpy record arrays ==

1. Database cells

Each databaes "cell" is a giant record array (on the order of a gigabyte for several million features SIFT, each 1024 bits long). Each "row" has the dtype

  sift_dtype = {
    'names': ['vec', 'geom', 'index'],
    'formats': ['128uint8', '4float32', 'uint16'],
  }

The 'vec' segment is the 1024 bit SIFT feature, and the 'geom' segment is the x,y,rotation, and hessian of feature.

Notice that the file from which the vector was extracted is not present! This is because storing filenames for each vector would be highly inefficient. Instead, they have a 16-bit index value, which indexes into another vector called the mapping. This is a map from indices->filenames.

2. Pixel to 3dspace maps

[pixels.py] proves a way to find the 3d location of a feature (as fetched from earthmine). While it caches lookups on disk, it does not provide high throughput when accessing random features in a cell, since it has to load a file for each different image the feature is from.

To get faster access to this data (for example, if you want to use the ratio test, but only test distance scores against feature scores spatially distant: see spatial-ratio-test in query.py), use the load_3dmap_for_cell method in the reader class. This returns a packed vector where the index of the feature in the cell vector corresponds directly to the index of the 3d point in this vector.

See the file for more specific usage instructions.

== The kdtree indexes ==

Another class of large files you will (maybe) encounter are the kdtree indexes built by FLANN. Be aware that if you want to conduct ANN searches outside of [query.py], you will have to manually save these indexes for fast searches.

== Configuration ==

[config.py] contains some deprecated options, you need not worry about them. However, the INFO(x) and INFO_TIMING(x) functions provide an easy way to output debug messages. They print time and thread information.

When running the query system, the following environment variables may apply:

	DEBUG: if DEBUG is set in the environment, only a single process will be used to do homography calcuations/draw images (to avoid intermingling of messages and lost exceptions). Note that to enable multiprocessing wrapping the characterize statement with a "with" statement:

		with system.Multiprocessing():
			system.characterize(C)
	
	Note that this is not related to the number of threads used when running database queries, which is managed by estimate_threads_avail() in [system.py].

	NO_HOM: if this is set the the query system will only compute the topN images matches and will not proceed with further processing.

	NO_DRAW: if this is set the query system will do homography/pose estimation logic, but will not actually draw an output imgae.

	See [context.py] for exact details and other run options.

== SIFT ==

The file:
[reader.py] Provides fast read/write access to the SIFT features. The sift features are extracted using some binary available on the Internet with the default settings. You must do this beforehand, which involves converting the image files to the appropriate size and format. Once you have extracted the SIFT features into the txt files, 

    get_reader(typehint) returns the correct reader class for the descriptors you want to use. These reader objects abstract away caching and packing of the extracted individual sift files into a numpy array. The first time you build a database cell using reader.load_cell(dir), you should get message from the reader "Reading features..." which will take a while.

    FeatureReader is returned by get_reader(), and has most of the methods for accessing database files. See [query.py] for how the reader is used. There are reader files for other descriptors besides sift. You can try these out by setting C.params['descriptor'] to 'chog' or 'surf', but these probably won't work out of the box. The performance with other descriptors is typically much lower, so the code paths are less maintained.

    PointToViewsMap provides efficient lookups from lat/lon to earthmine views. It is also used for tag occlusion detection, and is slightly more reliable that OcclusionSummary, though not completely accurate either. This might be because we only have 3d points for feature points, not all the points in the images.

== Data Types ==

[context.py] contains the Query and Context objects.

In general, Query objects are named Q and contain gps location, sensor data, image path, etc. Context objects C contain configuration information about what database we are using, the set of Q objects to iterate over, etc.

The Query object contains the pgm_scale parameter, which is used by the tag drawing code to determine how to scale the coordinate system to the output image. It should be the ratio of the size of .pgm file the sift features were extracted from and the jpg image the tags are to be drawn onto. This parameter is set by the Context object in a few of the query sets already, so if you aren't sure what it should be look at those examples.

Note that individual database images have no representation besides their sift feature-file name. (This is probably unfortunate). Most commonly we refer to them as "siftname" or "matchedimg". These are usually of the form [lat],[lon],[view].sift.txt

== Intermediate result reprentations ==

When you run the query system, it can produce a variety of different files as output. he most common is the .res file, which are put into the C.matchdir directory. The directory names are of the form "matchescells(g=100,r=d=236.6),queryX,{params}", so that each unique query will produce a distinct directory. This enables future runs with the same parameter to reuse the same cached files, saving lots of time.

The query will also output a a .res-detailed.npy file, which contains the feature matches as well as the count of feature votes. These are saved in a numpy packed file, which can be loaded like so:

	import numpy
	numpy.load("file.npy")

The format of these files:
	[[imagename1, [{'db': (x,y), 'query': (x2,y2)}, {'db'...}, {'db'...}, ...]]
	 [imagename2, [{'db': (x,y), 'query': (x2,y2)}, {'db'...}, {'db'...}, ...]]
	 [imagename3, [{'db': (x,y), 'query': (x2,y2)}, {'db'...}, {'db'...}, ...]]
	 ...]
	
The x,y coordinates specify the points where the features were matched. This file provides a lot more useful information that can be used in the geometric correspondence and tag transfer stages.

== Feature Correspondence Computation ==

[corr.py] has a bunch of random functions related to computing features matches and a first try at pose computation. I would not recommend working with these directly. Instead, use querySystem.py as an interface. You can use the C.selection parameter to operate on a single image or subset of images in a queryset.

Important functions:

    rematch(C, Q, dbsift) runs a linear NN search on Q and the db image. You can depend on this to be fairly fast, though not fast enough for hundreds of queries.

	getSpatiallyOrdered(...) imposes a 1 dimensional spatial ordering constraint on matches by finding the longest increasing subsequence via standard dynamic programming techniques. This hasn't given much performance improvement on top of other techniques so it's not used. I believe there is prior work about spatial constraints based on rotational ordering might work better (but is significantly more complex).

    find_corr(matches, hom=False, ....) handles all homography and fundamental matrix calculations. It can filter matches by rotation and attempts multiple parameters until a good relation is found.

    isHomographyGood(H) is a fairly reliable measure of if the computed homography is sane (not upside down, scaled badly, etc).

    draw_matches(C, Q, matches, rsc_matches, ...)

    The CameraModel and compute_pose is a try at greedy search to minimize reprojection error. It currently works, but is not very good at optimizing and only minimizes lat/lon error.

== Tags ==

[tags.py] contains several classes related to tags

    Tag holds all data about a tag and many nice utilities. For example, tag.isVisible(source) is one of the tag occlusion detection methods.

    TagCollection initializes a set of tags from a tags.csv and bearings.csv file. It matches up data from both files so that bearings.csv (which contains the normal vector of tags to their plane (tagged manually)) can be maintained seperately from the earthmine-viewer compatible tags.csv or left out altogether.

    OcclusionSummary holds earthmine information about a point. It can be queried to check what views can see a point. Unfortunately the occlusion data from earthmine is not reliable or available in all cases.

== Image Rendering ==

[render_tags.py] 

    ImageInfo is intended to be a class describing a database image. There are some functions that required you to create this class. It unifies various db image sources such as cell phones, earthmine truck, etc.

    TaggedImage: provides ways for tags to be project onto images. The basic tag drawing methods are contained here, as are all the occlusion detection / tag culling algorithms.

	To transfer tags onto an image, do the following.
	(NOTE that corr.py does these already!)

	1. Get the source of tags. By default, context.tags will is a TagCollection of berkeley tags with normal bearings.

	2. Call TaggedImage.map_tags_<method>, where tagged image is a database image. This will return tags associated with their xy coords in the database image. There are a few methods of transferring tags with different benefits/drawbacks documented in [render_tags.py] file.

	3. Use a method of projecting tags onto your image. For example, multiply each point by a homography matrix, or send the points to the mobile device to be draw, etc.

== Database Utils ==

[android.py] Provides a way to read ImageoTag xml metadata and images.

[earthMine.py] Earthmine API calls.
