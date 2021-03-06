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

config.hsv_enabled = False
C.resultsdir = os.path.expanduser('~/topmatches_big')

for dist, e, a, amb in [
    ('gauss', 1, (0,0), 25),
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
          C.ambiguity = amb
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
