import queryContext as context

context.QUERY = 'query4'
context.params.update({
  'checks': 1024,
  'trees': 1,
  'vote_method': 'filter',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})

context.num_images_to_print = 1
context.corrfilter_printed = 0
context.put_into_dirs = 0
context.showHom = 1
context.vars_init()
context.characterize()
