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
C.QUERY = 'q5-test'
#C.QUERY = 'query5horizontal'
C.params.update({
  'algorithm': 'kdtree',
  'checks': 1024,
  'trees': 1,
  'vote_method': 'filter',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})
C.resultsdir = os.path.expanduser('~/topmatches')

### shuffle cells experiment
#C.overlap_method = 0.0
#C.params['confstring'] = 'shuffle_cells_%s_overlap' % str(C.overlap_method)
#C.shuffle_cells = True
C.ncells = 8
C.ambiguity = 75
C.amb_cutoff = C.ambiguity
config.hsv_enabled = False
### shuffle cells experiment

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

#### weighted union experiment
#C.one_big_cell = True
##config.using_weighted_union = lambda: True
#config.hsv_enabled = False
##reader.config_mem_pin(True)
#C.ncells = 9
#C.ambiguity = 75
#### weighted union experiment

if 'DEBUG' in os.environ:
  system.characterize(C)
else:
  with system.MultiprocessExecution():
    system.characterize(C)

# vim: et sw=2
