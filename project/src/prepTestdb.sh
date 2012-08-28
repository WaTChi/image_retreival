#!/bin/bash

if [ -d "testdb" ]; then
	rm -rf testdb
fi

mkdir testdb
cp tutorial/example_db_images/* testdb
python prepTestdb.py testdb 1 png
mkdir testdb/cells-236.6
python prepTestdb.py testdb 2 0
