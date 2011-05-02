== Overview ==

We project tags from 3d space onto database images. Then, we move them onto query images using the computed homography.

== Data Types ==

[context.py] contains the Query and Context objects.

In general, Query objects are named Q and contain gps location, sensor data, image path, etc. Context objects C contain configuration information about what database we are using, the set of Q objects to iterate over, etc.

Note that individual database images have no representation besides their sift feature-file name. (This is probably unfortunate).

== Feature Correspondence Computation ==

[corr.py] has a bunch of random functions related to computing features matches and a first try at pose computation. I would not recommend working with these directly. Instead, use querySystem.py as an interface. You can use the C.selection parameter to operate on a single image or subset of images in a queryset.

Important functions:

    rematch(C, Q, dbsift) runs a linear NN search on Q and the db image.

    find_corr(matches, hom=False, ....) handles all homography and fundamental matrix calculations. It can filter matches by rotation and attempts multiple parameters until a good relation is found.

    isHomographyGood(H) is a fairly reliable measure of if the computed homography is sane (not upside down, scaled badly, etc).

    draw_matches(C, Q, matches, rsc_matches, ...)

    The CameraModel and compute_pose is a try at greedy search to minimize reprojection error. It currently works, but is not very good at optimizing and only minimizes lat/lon error.

== Tags ==

[tags.py] contains several classes related to tags

    Tag holds all data about a tag and many nice utilities.

    TagCollection initializes a set of tags from a tags.csv and bearings.csv file.

    OcclusionSummary holds earthmine information about a point. It can be queried to check what views can see a point. Unfortunately the occlusion data from earthmine is not reliable or available in all cases.

== Image Rendering ==

[render_tags.py] 

    ImageInfo is intended to be a class describing a database image. There are some functions that required you to create this class. It unifies various db image sources such as cell phones, earthmine truck, etc.

    TaggedImage: provides ways for tags to be project onto images. The basic tag drawing methods are contained here, as are all the occlusion detection / tag culling algorithms.

== Database Utils ==

[android.py] Provides a way to read ImageoTag xml metadata and images.

[reader.py] Provides fast read/write access to the SIFT features.

    get_reader(typehint) returns the correct reader class for the descriptors you want to use

    FeatureReader is returned by get_reader(), and has most of the methods for accessing database files. See [query.py] for how the reader is used.

    PointToViewsMap provides efficient lookups from lat/lon to earthmine views. It is slightly more reliable that OcclusionSummary.


