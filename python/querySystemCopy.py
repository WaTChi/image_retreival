import queryContext as context

context.QUERY = 'query3'
context.params.update({
  'checks': 1024,
  'trees': 1,
  'vote_method': 'filter',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})

context.vars_init()
context.characterize()
