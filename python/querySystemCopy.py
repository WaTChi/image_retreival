import system
from objects import DEFAULT_CONTEXT

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

#maindir = os.path.expanduser('~/.gvfs/data on gorgan.eecs.berkeley.edu/')
C.maindir = '/media/00C8173649E2CB4C/jz'

C.num_images_to_print = 1
C.corrfilter_printed = 0
C.do_posit = 1
C.put_into_dirs = 0
C.showHom = 1

system.characterize(C)
