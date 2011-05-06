#!/usr/bin/env python

print 'import modules...',
import sys
sys.stdout.flush()
print

import system
import os
from context import DEFAULT_CONTEXT

C = DEFAULT_CONTEXT.copy()
C.QUERY = 'query1'
C.params.update({
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
C.do_posit = 0
C.solve_pnp = 0
C.compute2dpose = 0 # [experimental, not recommended]
C.dump_hom = 0
C.ranking_min_consistent = 1

## Query3 tagging issues
#C.selection = ['8842', '8846', '8853', '8860', '8889', '8926']

# Query2 tagging issues
#C.selection = ['7727', '7735', '7744', '7746', '7751', '7753', '7755', '7756', '7763', '7764', '7765', '7776', '7710', '7712', '7716', '7717', '7739', '7747', '7764']

if 'DEBUG' in os.environ:
  system.characterize(C)
else:
  with system.MultiprocessExecution():
    system.characterize(C)

# vim: et sw=2
