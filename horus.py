import os
import subprocess

import util

__author__ = "zhangz"
__date__ = "$Dec 1, 2010 6:40:52 PM$"


def getFeatures(celldirpath):

    cells = util.getdirs(celldirpath)

    featureExtrator = "/home/zhangz/Downloads/libpmk-2.5/libpmk_features/tools/jz-images-to-pointsets.out"
    featureMapper = "/home/zhangz/Downloads/libpmk-2.5/libpmk_features/tools/jz-images-to-pointsets-mapping.out"
    for cell in cells:
        cellpath = os.path.join(celldirpath, cell)
        outputFilePath = os.path.join(celldirpath, cell + '.psl')
        outputMapFilePath = os.path.join(celldirpath, cell + '.map')
        if  not os.path.exists(outputFilePath):
            print "creating psl for cell: {0}".format(cell)
            subprocess.call([featureExtrator, outputFilePath, cellpath], cwd='/home/zhangz/Downloads/libpmk-2.5/libpmk_features/tools')
        if  not os.path.exists(outputMapFilePath):
            print "mapping psl for cell: {0}".format(cell)
            subprocess.call([featureMapper, outputMapFilePath, cellpath], cwd='/home/zhangz/Downloads/libpmk-2.5/libpmk_features/tools')

    cluster_builder = "/home/zhangz/Downloads/libpmk-2.5/libpmk2/tools/hierarchical-cluster-point-set.out"
    for cell in cells:
        pslpath = os.path.join(celldirpath, cell + '.psl')
        outputFilePath = os.path.join(celldirpath, cell + '.hc')
        if  not os.path.exists(outputFilePath):
            print "creating hc for cell: {0}".format(cell)
            subprocess.call([cluster_builder, pslpath, outputFilePath, str(4), str(10)], cwd='/home/zhangz/Downloads/libpmk-2.5/libpmk2/tools/')

    pyramid_builder = "/home/zhangz/Downloads/libpmk-2.5/libpmk2/tools/clusters-to-pyramids.out"
    for cell in cells:
        pslpath = os.path.join(celldirpath, cell + '.psl')
        hcpath = os.path.join(celldirpath, cell + '.hc')
        outputFilePath = os.path.join(celldirpath, cell + '.mrh')
        if  not os.path.exists(outputFilePath):
            print "creating mrh for cell: {0}".format(cell)
            subprocess.call([pyramid_builder, pslpath, hcpath, outputFilePath], cwd='/home/zhangz/Downloads/libpmk-2.5/libpmk2/tools/')

if __name__ == "__main__":
    getFeatures('/media/data/Research/cellsg=50,r=100,d=86.6')
