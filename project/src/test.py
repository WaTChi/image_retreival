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
import computePose
from numpy import transpose as tp

#if __name__ == '__main__':
#
#    imgpath = '/media/DATAPART2/ah/vp-test.jpg'
#    out0 = '/media/DATAPART2/ah/out3.jpg'
#    out1 = '/media/DATAPART2/ah/out1.jpg'
#    out2 = '/media/DATAPART2/ah/out2.jpg'
#    fov = 90 * np.pi/180
#    K = geom.cameramat(2500, 1200, fov)
#
#    answer = vanPts.vanishing_points(imgpath,K)
#
#    #vps, errs, numi = vanPts.vanishing_points(imgpath,K)
#
#    #print np.append( errs[:,np.newaxis] , np.append( numi[:,np.newaxis] , vps , 1 ) , 1 )


if __name__ == '__main__':

    normal = geom.normalrows(np.array([10*rnd.randn(),rnd.randn(),10*rnd.randn()]))
    offset = 10*abs(rnd.randn())

    x = tp([100*rnd.randn(100),100*rnd.randn(100),np.zeros(100),np.ones(100)])

    peq = np.append(normal,-offset)

    x[:,2] = -np.inner(x,peq) / peq[2]

    x[:,:3] += .1*offset * rnd.randn(100,3)

    n, pd = computePose.testplanefrom3d(x)

    print 'Plane equation: ' + str(peq)
    print 'Normal computed: ' + str(n)
    print 'Plane distance: %f' % pd
