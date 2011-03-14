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

C.max_matches_to_analyze = 5
C.stop_on_homTrue = 1
C.put_into_dirs = 1
C.do_posit = 0
C.dump_hom = 0

with system.MultiprocessExecution():
  system.characterize(C)

# vim: et sw=2
