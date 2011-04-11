#!/usr/bin/env python

import system
from context import DEFAULT_CONTEXT

C = DEFAULT_CONTEXT.copy()
C.QUERY = 'query2'
C.params.update({
  'checks': 1024,
  'trees': 1,
  'vote_method': 'filter',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})

C.max_matches_to_analyze = 1
C.stop_on_homTrue = 1
C.put_into_dirs = 0
C.do_posit = 0
C.solve_pnp = 0
C.compute2dpose = 0
C.dump_hom = 0
C.ransac_max_filt = 20
C.ransac_min_filt = 1

## Query3 tagging issues
#C.selection = ['8842', '8846', '8853', '8860', '8889', '8926']

# Query2 tagging issues
C.selection = ['7727', '7735', '7744', '7746', '7751', '7753', '7755', '7756', '7763', '7764', '7765', '7776']

with system.MultiprocessExecution():
  system.characterize(C)

# vim: et sw=2
