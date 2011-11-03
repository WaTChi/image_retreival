#!/usr/bin/env python

import os
import numpy as np
import reader
import sys
from config import *

# for running q5
# union of 15
big_cellpath = '/media/DATAPART2/Research/cells/g=100,r=d=236.6'
big_cell = '/media/DATAPART2/Research/cells/g=100,r=d=236.6/largeunion:15:116963924352de79e-sift.npy'
big_cell_pydata = '/media/DATAPART2/Research/cells/g=100,r=d=236.6/largeunion:15:116963924352de79e-sift-pydata.npy'
big_cell_map3d = '/media/DATAPART2/Research/cells/g=100,r=d=236.6/largeunion:15:116963924352de79e-map3d.npy'
pixmap_dir = '/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829'

print "Loading"
all_features = np.load(big_cell)
mapping = np.load(big_cell_pydata).item()
#reader.get_reader('sift').load_3dmap_for_cell(
#    big_cellpath,
#    all_features,
#    mapping,
#    pixmap_dir)

num_cells = 15
row_dtype = reader.SIFTReader.sift_dtype

print '%d features total' % len(all_features)
for overlap_percent in [0.0, 0.25, 0.5, 1.0, 1.50, 2.0, 3.0]:
    out_path = '/media/DATAPART2/Research/cells/test_artificial_overlap_%f' % overlap_percent
    if not os.path.isdir(out_path):
        os.mkdir(out_path) 
    cell_size = len(all_features)/15 * (1 + overlap_percent)
    print 'cell overlap %f' % (cell_size/float(len(all_features))*15)
    print 'cell size', cell_size

    for i in range(num_cells):
        celldir = os.path.join(out_path, 'sample_%d' % i)
        cellpath = os.path.join(out_path, 'sample_%d-sift.npy' % i)
        datapath = os.path.join(out_path, 'sample_%d-sift-pydata.npy' % i)
        map3d = os.path.join(out_path, 'sample_%d-map3d.npy' % i)
        if not os.path.exists(cellpath):
            print "Building", cellpath
            dataset = np.ndarray(cell_size, row_dtype)
            lookup_table = {}
            for i in range(len(dataset)):
                if i % 10000 == 0:
                    print '\r>> ' + str(i),
                    sys.stdout.flush()
                offset = int(i*len(all_features)/num_cells)
                dataset[i] = all_features[(offset + i) % len(all_features)]
            save_atomic(lambda d: np.save(d, dataset), cellpath)
        if not os.path.exists(datapath):
            os.symlink(big_cell_pydata, datapath)
#        if not os.path.exists(map3d):
#            os.symlink(big_cell_map3d, map3d)
        if not os.path.exists(celldir):
            os.mkdir(celldir)

