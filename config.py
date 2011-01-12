import os
import random
from datetime import datetime

LOCAL_CACHING = True
CACHE_PATH = os.path.expanduser('~/cache')
IS_REMOTE = lambda d: LOCAL_CACHING and '.gvfs' in d

def save_atomic(save_f, dest):
  if not os.path.exists(CACHE_PATH):
      os.mkdir(CACHE_PATH)
  a = os.path.dirname(dest)
  b = os.path.basename(dest)
  tmp = os.path.join(a, '~%d' % random.randint(1e4,1e5) + b)
  try:
    save_f(tmp)
    os.rename(tmp, dest)
  finally:
    if os.path.exists(tmp):
      os.unlink(tmp)

def INFO(x):
  print datetime.today().strftime("%l:%M:%S - ") + str(x)

def INFO_TIMING(x):
#  INFO(x)
  pass

# vim: et sw=2
