__author__="aaronh"
__date__ ="$Jul 25, 2011 12:46:02 PM$"

from distutils.core import setup, Extension

setup (
  name = 'fib',
  version = '1.0',

  ext_modules=[Extension('fib', ['lsd.c'])],

  # Fill in these to make your Egg ready for upload to
  # PyPI
  author = 'aaronh',
  author_email = '',

  url = '',
  license = '',
  long_description= 'Long description of the package',

  # could also include long_description, download_url, classifiers, etc.

  
)
