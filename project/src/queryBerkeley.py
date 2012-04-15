#!/usr/bin/env python

print 'import modules...',
import sys
sys.stdout.flush()
print

import system
import reader
import os
from context import DEFAULT_CONTEXT

C = DEFAULT_CONTEXT.copy()
C.QUERY = 'query5horizontal'
C.params.update({
  'algorithm': 'kdtree',
  'checks': 1024,
  'trees': 1,
  'vote_method': 'filter',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})

C.max_matches_to_analyze = 1
C.stop_on_homTrue = 0
C.put_into_dirs = 0
C.show_feature_pairs = True
C.do_posit = 0
C.solve_pnp = 0
C.dump_hom = 0

C.ranking_min_consistent = 1
C.ranking_max_considered = 30

C.ncells = 9
C.ambiguity = 75

# Pose estimation parameters

C.hiresdir='/media/DATAPART2/Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829/hires'
C.solve_pose = 1
C.pose_remove = [ '2011-04-04_15-50-21_519' , '2011-04-04_15-22-31_544' ,
                  '2011-04-04_14-56-08_926' , '2011-04-04_15-06-58_888' ,
                  '2011-04-04_15-19-11_784' , '2011-04-04_15-32-55_864' ,
                  '2011-04-04_15-36-15_712' , '2011-04-04_15-06-48_388' ,
                  '2011-04-04_15-23-24_958' , '2011-04-04_14-57-13_794' ,
                  ### QUERIES REMOVED BECAUSE OF MISSING DEPTH IMAGE ###
                  '2011-04-04_15-02-39_082' , '2011-04-04_15-06-37_147' , 
                  '2011-04-04_15-10-56_749' , '2011-04-04_15-15-08_854' ,
                  '2011-04-04_15-17-45_763' , '2011-04-04_15-19-27_080' , 
                  '2011-04-04_15-19-41_314' , '2011-04-04_15-35-13_268' ,
                  '2011-04-04_15-36-56_430' ]
#                  # STRICT ELIMINATION #
#C.pose_remove = [ '2011-04-04_15-50-21_519' , '2011-04-04_15-22-31_544' ,
#                  '2011-04-04_14-56-08_926' , '2011-04-04_15-06-58_888' ,
#                  '2011-04-04_15-19-11_784' , '2011-04-04_15-32-55_864' ,
#                  '2011-04-04_15-36-15_712' , '2011-04-04_15-35-13_268' ,
#                  '2011-04-04_15-23-24_958' , '2011-04-04_15-06-48_388' ,
#                  '2011-04-04_15-34-01_368' , '2011-04-04_15-47-38_651' ,
#                  '2011-04-04_15-19-41_314' , '2011-04-04_15-15-55_363' ,
#                  ### QUERIES REMOVED BECAUSE OF MISSING DEPTH IMAGE ###
#                  '2011-04-04_15-02-39_082' , '2011-04-04_15-06-37_147' , 
#                  '2011-04-04_15-10-56_749' , '2011-04-04_15-15-08_854' ,
#                  '2011-04-04_15-17-45_763' , '2011-04-04_15-19-27_080' , 
#                  '2011-04-04_15-19-41_314' , '2011-04-04_15-35-13_268' ,
#                  '2011-04-04_15-36-56_430' ]
#                  # STRICT ELIMINATION #
C.pose_param = {
    'resultsdir'    : '/media/DATAPART2/ah/pose_runs/berkeley',
    'run_info'      : '/media/DATAPART2/ah/pose_runs/berkeley/run_info.txt',
    'pose_file'     : '/media/DATAPART2/ah/pose_runs/berkeley/pose_results.txt',
    'extras_file'   : '/media/DATAPART2/ah/pose_runs/berkeley/extras.txt',
    'cheat'         : False,
    'runflag'       : 11,
    'remove_ground' : False,
    'use_weight'    : True,
    'solve_bad'     : True,
    'maxmatch'      : 1,
    'maxratio'      : 0.8,
    'maxdist'       : 10**7,
    'ransac_iter'   : 10**7,
    'inlier_error'  : 0.01 }


# SUBSELECTIONS
#C.selection = ['2011-04-04_15-32-28_500']
#C.selection = ['2011-04-04_15-20-49_89']

# Run in debug mode or parallel
debug = True
try: print 'Number of query images selected = %d' % len(C.selection)
except TypeError: debug=False
if 'DEBUG' in os.environ or debug:
  system.characterize(C)
else:
  with system.MultiprocessExecution():
    system.characterize(C)

# vim: et sw=2
