import queryContext as context

def run():
    context.QUERY = 'query4-matlab'
    context.params.update({
      'checks': 1024,
      'trees': 1,
      'vote_method': 'filter',
      'num_neighbors': 1,
      'dist_threshold': 70000,
      'confstring': '',
    })
    context.print_per = 1000
    context.ambiguity = 75
    context.topnresults = [1,2,5,10]
#    context.locator_function = context.skew_location
    context.locator_function = context.load_location
    context.match_callback = context.dump_combined_matches
    context.cacheEnable = 1
    context.ransac_min_filt = 100
    context.ncells = 10
    context.vars_init()
    context.characterize()

    for n in context.topnresults:
        total_count = context.results[n]
        match_rate = float(total_count) / context.count
        print "matched {0}\t out of {1}\t = {2}\t in the top {3}".format(total_count, context.count, match_rate, n)

run()
