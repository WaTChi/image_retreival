#!/usr/bin/env python

import system
from context import DEFAULT_CONTEXT

C = DEFAULT_CONTEXT.copy()
C.QUERY = 'query3'
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
C.dump_hom = 0
C.ransac_min_filt = 1
C.selection = ['8852']

with system.MultiprocessExecution():
  system.characterize(C)

# vim: et sw=2
