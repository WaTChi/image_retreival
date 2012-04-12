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
C.QUERY = 'oakland1'
C.params.update({
  'algorithm': 'kdtree',
  'checks': 1024,
  'trees': 1,
  'vote_method': 'top_n',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})

# cell compression criteria
#C.criteria = 'random_quarter'
#C.criteria = 'locally_significant'
#C.criteria = 'locally_insignificant'

config.hsv_enabled = False
#reader.config_mem_pin(True)
C.resultsdir = os.path.expanduser('~/topmatches_big')
C.one_big_cell = 0

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

for dist, e, a, co in [
    ('gauss', 1, (0,0), 75),
#    ('uniform', 0.5, (0, 50**2), 75),
#    ('uniform', 0.5, (0, 125**2), 150),
#    ('uniform', 0.5, (0, 225**2), 250),
#    ('uniform', 0.5, (0, 325**2), 350),
  ]:
  for seed in [0]:
    for r,d in [
#        (500, 500),
#        (350, 350),
        (236.6, 236.6),
#        (150, 150),
        ]:
      C._test_r = r
      C._test_d = d
      C.cellradius = r
      for nn, ncells in [(1,999)]:
          C.amb_cutoff = co
          C.ambiguity = co
          C.added_error = {
            'seed': seed,
            'dist': dist,
            'exponent': e,
            'args': a,
            'amb_padding': 25,
          }
          C.amb_padding = 25
          C.ncells = ncells
          C.params['num_neighbors'] = nn
          C.params['trees'] = 1
          with system.MultiprocessExecution():
            system.characterize(C)

# vim: et sw=2
