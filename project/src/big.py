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
  'checks': 1024,
  'trees': 4,
  'vote_method': 'top_n',
  'num_neighbors': 4,
  'dist_threshold': 70000,
  'confstring': '',
})
config.hsv_enabled = False
#reader.config_mem_pin(True)
C.resultsdir = os.path.expanduser('~/topmatches_big')
C.ambiguity = 75

# q5h std=9.5 (var=91.8), so to get std=25m (var=625)
# we add 23 (sqrt(625-91.8)) as variance
C.added_error = {
  'seed': 0,
  'dist': 'gauss',
  'args': (0, 23),
}
C.amb_cutoff = C.ambiguity
C.ncells = 4
C.one_big_cell = 0

#C.restrict_cells = set(['37.8732916024,-122.265650441', '37.8696062215,-122.273737308', '37.8732916946,-122.273737207', '37.871448912,-122.269693992', '37.8677634389,-122.269694228', '37.8696061908,-122.271041854', '37.8696061293,-122.265650946', '37.8732916639,-122.271041618', '37.8714489427,-122.272389514', '37.8714488505,-122.26430295', '37.8732916331,-122.268346029', '37.8714489734,-122.275085035', '37.8714488812,-122.266998471', '37.8677634696,-122.272389615', '37.86960616,-122.2683464'])

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

for seed in range(3):
  for nn in [1,2,4,8,16]:
    for ncells in [1,2,4,8]:
      C.added_error = {
        'seed': seed,
        'dist': 'gauss',
        'args': (0, 23),
#        'dist': 'uniform',
#        'args': (0, 50),
      }
      C.ncells = ncells
      C.params['num_neighbors'] = nn
      with system.MultiprocessExecution():
        system.characterize(C)

# vim: et sw=2
