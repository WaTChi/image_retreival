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
C.QUERY = 'oakland1'
C.params.update({
  'algorithm': 'kdtree',
  'checks': 1024,
  'trees': 1,
  'vote_method': 'top_n',
  'num_neighbors': 8,
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

C.ncells = 1
C.ambiguity = 75

# Pose estimation parameters

C.hiresdir='/media/DATAPART1/oakland/earthmine/rect_hires'
C.solve_pose = 1
C.pose_remove = [ '2011-10-28_12-45-09_775' , '2011-10-28_12-23-24_368' ,
                  '2011-10-28_12-13-41_589' , '2011-10-28_12-12-29_132' ,
                  '2011-10-28_12-12-20_282' , '2011-10-28_12-11-42_811' ,
                  '2011-10-28_12-09-45_828' , '2011-10-28_12-09-35_578' ,
                  '2011-10-28_12-01-35_821' , '2011-10-28_11-59-53_593' ]
#                  # STRICT ELIMINATION #
#C.pose_remove = [ '2011-10-28_12-45-09_775' , '2011-10-28_12-23-24_368' ,
#                  '2011-10-28_12-13-41_589' , '2011-10-28_12-12-29_132' ,
#                  '2011-10-28_12-12-20_282' , '2011-10-28_12-11-42_811' ,
#                  '2011-10-28_12-09-45_828' , '2011-10-28_12-09-35_578' ,
#                  '2011-10-28_12-01-35_821' , '2011-10-28_11-59-53_593' ,
#                  '2011-10-28_13-00-41_047' , '2011-10-28_12-50-32_083' ,
#                  '2011-10-28_12-01-26_932' , '2011-10-28_11-59-28_617' ]
#                  # STRICT ELIMINATION #
C.pose_param = {
    'resultsdir'    : '/media/DATAPART2/ah/pose_runs/oakland',
    'run_info'      : '/media/DATAPART2/ah/pose_runs/oakland/run_info.txt',
    'pose_file'     : '/media/DATAPART2/ah/pose_runs/oakland/pose_results.txt',
    'extras_file'   : '/media/DATAPART2/ah/pose_runs/oakland/extras.txt',
    'cheat'         : False,
    'runflag'       : 11,
    'remove_ground' : True,
    'use_weight'    : True,
    'solve_bad'     : False,
    'maxmatch'      : 1,
    'maxratio'      : 0.8,
    'maxdist'       : 10**7,
    'ransac_iter'   : 10**7,
    'inlier_error'  : 0.01 }


# SUBSELECTIONS
#C.selection = ['2011-10-28_13-00-41_047']
#C.selection = ['2011-10-28_11-59-08_593']
#C.selection = ['2011-10-28_11-59-08_593']
### For checking dyaw weighting ###
#C.selection = ['2011-10-28_12-56-32_511','2011-10-28_12-39-00_773','2011-10-28_12-02-44_400','2011-10-28_12-02-23_300']

debug = True
try: print 'Query selection size: %d' % len(C.selection)
except TypeError: debug = False
if 'DEBUG' in os.environ or debug:
  system.characterize(C)
else:
  with system.MultiprocessExecution():
    system.characterize(C)

# vim: et sw=2
