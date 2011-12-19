#!/usr/bin/env python

print 'import modules...',
import config
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
  'checks': 4096,
  'trees': 1,
  'vote_method': 'top_n',
  'num_neighbors': 4,
  'dist_threshold': 70000,
  'confstring': '',
})
C.resultsdir = os.path.expanduser('~/topmatches_x8')

config.hsv_enabled = False
C.ncells = 1
C.ambiguity = 75

C.max_matches_to_analyze = 1
C.stop_on_homTrue = 0
C.put_into_dirs = 0
C.show_feature_pairs = True
C.do_posit = 0
C.solve_pnp = 0
C.solve_pose = 0
C.solve_bad = 0
C.compute2dpose = 0 # [experimental, not recommended]
C.dump_hom = 0

C.spatial_comb = 0 # doesn't work very well

C.ranking_min_consistent = 1
C.ranking_max_considered = 30
C.weight_by_distance = False
C.weight_by_coverage = False

if 'DEBUG' in os.environ:
  system.characterize(C)
else:
  with system.MultiprocessExecution():
    system.characterize(C)

# vim: et sw=2
