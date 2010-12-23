import os

LOCAL_CACHING = True
CACHE_PATH = os.path.expanduser('~/cache')
IS_REMOTE = lambda d: LOCAL_CACHING and '.gvfs' in d

def INFO(x):
  print 'INFO: ' + str(x)
