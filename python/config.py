import sys
import os
import random
import threading
import shutil
from datetime import datetime

LOCAL_CACHING = False
CACHE_PATH = os.path.expanduser('/var/tmp')
IS_REMOTE = lambda d: LOCAL_CACHING and '.gvfs' in d or '/shiraz/' in d
DETAIL_VERSION = str(15)

if os.getenv('SERVER_PROTOCOL'):
  STDOUT = sys.stderr
else:
  STDOUT = sys.stdout

def getdests(directory, name):
  if type(directory) is list:
    directory = directory[0]
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
  if type(directory) is list:
    directory = directory[0]
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

def getcellid(cellpath):
  if type(cellpath) is list:
    cellpath.sort()
    cellid = 'union:' + ':'.join([os.path.basename(d) for d in cellpath])
    if len(cellid) > 200:
      cellid = 'largeunion:%d:%x' % (len(cellpath), abs(sum(map(hash, cellpath))))
    return cellid
  return os.path.basename(cellpath)

def INFO(x):
  print >> STDOUT, "[%d]" % os.getpid(), datetime.today().strftime("%l:%M:%S - ") + str(x)

def INFO_TIMING(x):
#  INFO(x)
  pass

# vim: et sw=2
