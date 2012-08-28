#!/bin/bash

rm -rf corydb
mkdir corydb
cp tutorial/cory_nikon/*.jpg corydb

python prepdb.py corydb 1 jpg 
python prepdb.py corydb 3 jpg

mkdir corydb/cells-100
python prepdb.py corydb 2 0

python prepdb.py tutorial/cory_nexus 1 jpg

