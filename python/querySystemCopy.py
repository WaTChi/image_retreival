#!/usr/bin/env python

import system
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

C.max_matches_to_analyze = 5
C.corrfilter_printed = 1
C.do_posit = 0
C.put_into_dirs = 1
C.dumphom = 0

# enable multiprocessing
C.pool_enable()
system.characterize(C)
C.pool_shutdown()
