This directory contains ALL Python code for image retrieval and pose estimation. While it should be cleaned up into subdirectories, no one has yet had the time to do that, so this documentation splits the functions into categories.

*** See Documentation Folder for details ***


QUERY WRAPPER FILES: These files call the image retrieval and pose estimation system. They look very similar. The procedure for running a new experiment is typically to copy one of the files and modify it.
	queryBerkeley
	queryBerkeleyCheat
	queryOakland
	queryOaklandCheat
    queryOakland
    querySystem
    querySystemCopy
    querySystemCory

	
IMAGE RETRIEVAL:
    See Documentation folder for details.

POSE ESTIMATION:
	computePose : 
		top-level function which organizes metadata and computes the pose for the query.
	solveEssMatrix :
		solves the constrained essential matrix (not used in pose estimation but still good to keep).
	solveHomography :
		solves the contrained homography matrix with varying knowns and unknowns (used in pose estimation).
	vp_analysis :
		performs the vanishing point detection and alignment for pose estimation.
	

UTILITY:
    See Documentation folder for details.

	extract_oakland:
		contains code which sets up the oakland database, including converting panoramas to rectilinear images.
	geom : 
		contains most of the geometry functions, like conversions of lat/lon to grid or yaw, pitch, roll to rotation matrices.
	util : 
		contains general utility functions for image retrieval.
		
		
GROUND TRUTH FILES: 
    Oakland:
        oakland1GroundTruth :
    Berkeley:
        groundtruth*.py
        query*GroundTruth.py

