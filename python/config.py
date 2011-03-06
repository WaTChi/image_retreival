import sys
import os
import random
import shutil
from datetime import datetime

LOCAL_CACHING = True
CACHE_PATH = os.path.expanduser('/var/tmp')
IS_REMOTE = lambda d: LOCAL_CACHING and '.gvfs' in d or '/shiraz/' in d

if os.getenv('SERVER_PROTOCOL'):
  STDOUT = sys.stderr
else:
  STDOUT = sys.stdout

def getdests(directory, name):
  default = os.path.join(os.path.dirname(directory), name)
  if IS_REMOTE(directory):
    INFO('copying %s to local cache' % name)
    local = os.path.join(CACHE_PATH, name)
    return [default, local]
  return [default]

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

def getfile(directory, name):
  default = os.path.join(os.path.dirname(directory), name)
  if IS_REMOTE(directory):
    local = os.path.join(CACHE_PATH, name)
    if os.path.exists(local):
      return local
    elif os.path.exists(default):
      INFO('copying %s to local cache' % name)
      save_atomic(lambda d: shutil.copyfile(default, d), local)
      return local
  return default

def getcellid(directory):
  return os.path.basename(directory)

def INFO(x):
  print >> STDOUT, datetime.today().strftime("%l:%M:%S - ") + str(x)

def INFO_TIMING(x):
#  INFO(x)
  pass

# vim: et sw=2
