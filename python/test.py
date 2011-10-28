#import computePose
import time
import reader
import pyflann
import numpy.linalg as alg
import numpy as np
import cv
import Image
import os
import numpy.random as rnd
import vanPts
import geom

if __name__ == '__main__':
    
    imgpath = '/media/DATAPART2/ah/vp-test.jpg'
    out0 = '/media/DATAPART2/ah/out3.jpg'
    out1 = '/media/DATAPART2/ah/out1.jpg'
    out2 = '/media/DATAPART2/ah/out2.jpg'
    fov = 90 * np.pi/180
    K = geom.cameramat(2500, 1200, fov)
    
    answer = vanPts.vanishing_points(imgpath,K)

    #vps, errs, numi = vanPts.vanishing_points(imgpath,K)

    #print np.append( errs[:,np.newaxis] , np.append( numi[:,np.newaxis] , vps , 1 ) , 1 )
    


