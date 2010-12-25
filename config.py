import os
from datetime import datetime

LOCAL_CACHING = True
CACHE_PATH = os.path.expanduser('~/cache')
IS_REMOTE = lambda d: LOCAL_CACHING and '.gvfs' in d

if not os.path.exists(CACHE_PATH):
    os.mkdir(CACHE_PATH)

def INFO(x):
  print datetime.today().strftime("%l:%M:%S - ") + str(x)

def INFO_TIMING(x):
#  INFO(x)
  pass

# vim: et sw=2
