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
C.solve_bad = 0
C.solve_pose = 1
C.pose_param = {
    'resultsdir'    : '/media/DATAPART2/ah/pose_runs/berkeley',
    'pose_file'     : '/media/DATAPART2/ah/pose_runs/berkeley/pose_results.txt',
    'mode'          : 'hom',
    'runflag'       : 6,
    'use_planes'    : False,
    'remove_ground' : False,
    'min_domsize'   : 0.5,
    'maxratio'      : 0.5,
    'maxdist'       : 50000,
    'ransac_iter'   : 1000,
    'ransac_depth'  : True,
    'inlier_error'  : 0.05 }


# SUBSELECTIONS
#C.selection = ['2011-04-04_15-42-43_581']
#C.selection = ['2011-04-04_15-32-55_864']

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
