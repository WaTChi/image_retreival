# Thread-safe allocator of flann instances and other large memory objects.
# Use this class for efficient repeated loads of large files.
# 
# Usage examples:
#
# from memory import Shared
# Shared.get_object(...)
# Shared.get_array(...)
# Shared.init_flann(...)

import numpy as np

class MemoryManager(object):
  __instance = None
  def __init__(self):
    if MemoryManager.__instance:
      raise Exception("Use the memory.Shared instance")
    MemoryManager.__instance = self

  def get_object(self, npyfile):
    pass

  def get_array(self, npyfile):
    pass

  def init_flann(self, celldir):
    pass

Shared = MemoryManager()

# vim: et sw=2
