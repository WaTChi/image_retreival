This directory contains the functions which perform analysis on the results of pose estimation. The analysis is meaningless without a correct ground truth, however.


FUNCTIONS:

pose_analysis.m:
	This function analyzes latitude and longitude error of the homography estimate with respect to the ground truth. Plots a performance curve as output to the function.
	(INPUTS)
	runnum - Which run to analyze (DEFAULT = 0)
	setnum - Which set to analyze, 5 = Berkeley, 11 = Oakland (DEFAULT = 11)
	gps	   - Boolean which plots the GPS performance with respect to ground truth if set to True (DEFAULT = 0)
	(EXAMPLES)
	pose_analysis -> Runs analysis on the results in /oakland in the pose estimation results directory
	pose_analysis(2,5) -> Runs analysis on the results in /berkeley2 in the pose estimation results directory.
	pose_analysis(42,11,1) -> Runs analysis on the results in /oakland42 and includes a plot of GPS performance.
	
prehom_analysis.m:
	This function analyzes yaw orientation error of the vanishing point alignment algorithm with respect to the ground truth. The output is again a performance plot, and by default a plot containing the performance of the cell phone compass is included.
	(INPUTS)
	runnum - Which run to analyze (DEFAULT = 0)
	setnum - Which set to analyze, 5 = Berkeley, 11 = Oakland (DEFAULT = 11)
	(EXAMPLES)
	pose_analysis(5) -> Runs analysis on the results in /oakland6 in the pose estimation results directory
	pose_analysis(0,5) -> Runs analysis on the results in /berkeley in the pose estimation results directory.