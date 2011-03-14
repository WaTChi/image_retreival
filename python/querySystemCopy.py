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

C.num_images_to_print = 1
C.corrfilter_printed = 0
C.do_posit = 1
C.put_into_dirs = 0
C.showHom = 1

system.characterize(C)
