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

C.hiresdir='/media/DATAPART1/oakland/earthmine/rect_hires'
C.solve_bad = 0
C.solve_pose = 1
C.pose_param = {
    'resultsdir'    : '/media/DATAPART2/ah/pose_runs/oakland',
    'pose_file'     : '/media/DATAPART2/ah/pose_runs/oakland/pose_results.txt',
    'mode'          : 'comb',
    'runflag'       : 6,
    'use_planes'    : True,
    'remove_ground' : True,
    'min_domsize'   : 0.3,
    'maxratio'      : 0.5,
    'maxdist'       : 50000,
    'ransac_iter'   : 1000,
    'ransac_depth'  : True,
    'inlier_error'  : 0.03 }


# SUBSELECTIONS
#C.selection = ['2011-10-28_11-51-29_558','2011-10-28_11-54-24_330','2011-10-28_11-55-43_140','2011-10-28_11-56-05_366','2011-10-28_11-56-53_831'] # query orientation check
#C.selection = ['2011-10-28_11-56-05_366']
#C.selection = ['2011-10-28_11-51-29_558'] # first image, check orientation
#C.selection = ['2011-10-28_12-46-56_127'] # 0 inliers, good pose?
#C.selection = ['2011-10-28_12-04-57_265'] # 0 inliers with result?
#C.selection = ['2011-10-28_11-57-59_152']
#C.selection = ['2011-10-28_13-00-12_019'] # doesn't work for some reason
C.selection = ['2011-10-28_11-55-43_140']

if 'DEBUG' in os.environ:
  system.characterize(C)
else:
  with system.MultiprocessExecution():
    system.characterize(C)

# vim: et sw=2
