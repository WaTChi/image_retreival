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
C.params.update({
  'algorithm': 'kdtree',
  'checks': 1024,
  'trees': 1,
  'vote_method': 'filter',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})
C.resultsdir = os.path.expanduser('~/topmatches_193')
#C.added_error = 250

#C.restrict_cells = set(['37.8732916024,-122.265650441', '37.8696062215,-122.273737308', '37.8732916946,-122.273737207', '37.871448912,-122.269693992', '37.8677634389,-122.269694228', '37.8696061908,-122.271041854', '37.8696061293,-122.265650946', '37.8732916639,-122.271041618', '37.8714489427,-122.272389514', '37.8714488505,-122.26430295', '37.8732916331,-122.268346029', '37.8714489734,-122.275085035', '37.8714488812,-122.266998471', '37.8677634696,-122.272389615', '37.86960616,-122.2683464'])
config.hsv_enabled = False
C.ncells = 8
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
