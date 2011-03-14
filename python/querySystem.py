#!/usr/bin/env python

import system
from context import DEFAULT_CONTEXT

C = DEFAULT_CONTEXT.copy()

def run():
    C.QUERY = 'query4-matlab'
    C.params.update({
      'checks': 1024,
      'trees': 1,
      'vote_method': 'filter',
      'num_neighbors': 1,
      'dist_threshold': 70000,
      'confstring': '',
    })
    C.print_per = 1000
    C.ambiguity = 75
    C.topnresults = [1,2,5,10]
#    C.locator_function = system.skew_location
    C.locator_function = system.load_location
    C.match_callback = system.dump_combined_matches
    C.cacheEnable = 1
    C.ransac_min_filt = 100
    C.ncells = 10
    results, count = system.characterize(C)

    for n in C.topnresults:
        total_count = results[n]
        match_rate = float(total_count) / count
        print "matched {0}\t out of {1}\t = {2}\t in the top {3}".format(total_count, count, match_rate, n)

run()
