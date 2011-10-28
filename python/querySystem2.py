#!/usr/bin/env python

import os
import system
from context import DEFAULT_CONTEXT

C = DEFAULT_CONTEXT.copy()

def output(C, results, count):
    for n in C.topnresults:
        total_count = results[n]
        match_rate = float(total_count) / count
        print "matched {0}\t out of {1}\t = {2}\t in the top {3}".format(total_count, count, match_rate, n)

def run():
    C.QUERY = 'query5vertical'
    C.params.update({
      'checks': 1024,
      'trees': 1,
      'vote_method': 'filter',
      'num_neighbors': 1,
      'dist_threshold': 70000,
      'confstring': '',
    })
    C.print_per = 1000
    C.verbosity = 0
    C.ambiguity = 75
    C.topnresults = [1,2,3,4,5,6,7,8,9,10]
#    C.locator_function = system.skew_location
#    C.locator_function = system.load_location
    C.aarondir='stopearly10m'
    C.match_callback = system.dump_combined_matches
    C.cacheEnable = True
    C.tagcompute = True
    C.show_feature_pairs = False
    C.ranking_min_consistent = max(C.topnresults)
    C.ranking_max_considered = 100
    C.ncells = 10
    C.resultsdir="/media/DATAPART2/jz/tmp"
    if 'DEBUG' in os.environ:
        results, count = system.characterize(C)
    else:
        with system.MultiprocessExecution():
            results,count = system.characterize(C)
    output(C, results, count)

#    for q in ['query1', 'query2', 'query3','query4','query5horizontal','query5vertical']:
#        C.QUERY = q
#        C.aarondir='stopearly10m'

#        C.aarondir='stopearly'
#        C.ranking_min_consistent = max(C.topnresults)
#        print "query: {0}, method: stopearly".format(q)
#        results, count = system.characterize(C)
#        output(C, results, count)

#        C.aarondir='maxhundred'
#        C.ranking_min_consistent = float('inf')
#        print "query: {0}, method: max100".format(q)
#        results, count = system.characterize(C)
#        output(C, results, count)




run()
