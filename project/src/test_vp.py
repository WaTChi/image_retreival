# test vanishing point analysis

import solveHomography
import computePose
import numpy as np
import numpy.linalg as alg
import numpy.random as rnd
from numpy import transpose as tp
import geom
import computePose2
import vp_analysis

if __name__ == '__main__':

    imgpath2 = '/media/DATAPART1/oakland/earthmine/rect_hires/37.794703,-122.257721-0002.jpg'
    lsdpath2 = '/media/DATAPART1/oakland/earthmine/rect_hires/lsd/37.794703,-122.257721-0002.lsd'

    imgpath3 = '/media/DATAPART1/oakland/earthmine/rect_hires/37.794703,-122.257721-0003.jpg'
    lsdpath3 = '/media/DATAPART1/oakland/earthmine/rect_hires/lsd/37.794703,-122.257721-0003.lsd'

    imgpath4 = '/media/DATAPART1/oakland/earthmine/rect_hires/37.794703,-122.257721-0004.jpg'
    lsdpath4 = '/media/DATAPART1/oakland/earthmine/rect_hires/lsd/37.794703,-122.257721-0004.lsd'

    Kcal = geom.cameramat(2500, 1200, 90*np.pi/180)

    midpts2, lineqs2, lengths2 = vp_analysis.LfromLSD(lsdpath2, imgpath2, Kcal)
    midpts3, lineqs3, lengths3 = vp_analysis.LfromLSD(lsdpath3, imgpath3, Kcal)
    midpts4, lineqs4, lengths4 = vp_analysis.LfromLSD(lsdpath4, imgpath4, Kcal)

    print midpts2
    print lineqs2
    print lengths2

    Rot2 = geom.RfromYPR(60,0,0)
    Rot3 = geom.RfromYPR(90,0,0)
    Rot4 = geom.RfromYPR(120,0,0)

    vp_threshold = 5

    vps2, conf2 = vp_analysis.VPfromRANSAC(midpts2, lineqs2, lengths2, Rot2, vp_threshold)
    vps3, conf3 = vp_analysis.VPfromRANSAC(midpts3, lineqs3, lengths3, Rot3, vp_threshold)
    vps4, conf4 = vp_analysis.VPfromRANSAC(midpts4, lineqs4, lengths4, Rot4, vp_threshold)

    wvps2 = tp( np.dot(Rot2,tp(vps2)) )
    wvps3 = tp( np.dot(Rot3,tp(vps3)) )
    wvps4 = tp( np.dot(Rot4,tp(vps4)) )

    print 'Image vanishing points:'
    print np.concatenate( (vps2,tp([conf2])) , 1 )
    print np.concatenate( (vps3,tp([conf3])) , 1 )
    print np.concatenate( (vps4,tp([conf4])) , 1 )

    print 'World vanishing points:'
    print np.concatenate( (wvps2,tp([conf2])) , 1 )
    print np.concatenate( (wvps3,tp([conf3])) , 1 )
    print np.concatenate( (wvps4,tp([conf4])) , 1 )

#    tvp = np.array([[-1,0,1],[1,0,1],[1,1,1],[0,1,1],[-100,0,1]])
#    ninl = np.array([100,150,50,150,200])
#    nout = 350
#    print 'Vanishing points: '
#    print tvp
#
#    lineqs = np.zeros((np.sum(ninl)+nout,3))
#    midpts = np.zeros((np.sum(ninl)+nout,3))
#    idx = 0
#    for i in range(len(tvp)):
#        pix = np.concatenate( (rnd.rand(ninl[i],2)-0.5,tp([np.ones(ninl[i])])) , 1 )
#        leq = np.zeros((ninl[i],3))
#        leq[:,0] , leq[:,1] = tvp[i,1]-pix[:,1] , pix[:,0]-tvp[i,0]
#        leq[:,2] = -np.sum(leq*pix,1)
#        leq = geom.normalrows(leq)
#        leq[:,2] += 0 * (rnd.rand(ninl[i])-0.5)
#        lineqs[idx:idx+ninl[i],:] = leq
#        midpts[idx:idx+ninl[i],:] = pix
#        idx += ninl[i]
#    pix1 = np.concatenate( (rnd.rand(nout,2)-0.5,tp([np.ones(nout)])) , 1 )
#    pix2 = np.concatenate( (rnd.rand(nout,2)-0.5,tp([np.ones(nout)])) , 1 )
#    leq = np.zeros((nout,3))
#    leq[:,0] , leq[:,1] = pix1[:,1]-pix2[:,1] , pix2[:,0]-pix1[:,0]
#    leq[:,2] = -np.sum(leq*pix1,1)
#    leq = geom.normalrows(leq)
#    leq[:,2] += 0.1 * (rnd.rand(nout)-0.5)
#    lineqs[idx:,:] = leq
#    midpts[idx:,:] = (pix1+pix2)/2.0
#    lengths = 0.01 + 0.1*rnd.randn(np.sum(ninl)+nout)**2
#
#    Rot = geom.RfromYPR(0,0,0)
#    vp_threshold = 5
#
##    i = 4
##    a,b = 450,650
##    lineqs, lengths = lineqs[a:b,:], lengths[a:b]
#    vps, conf = vp_analysis.VPfromRANSAC(midpts, lineqs, lengths, Rot, vp_threshold)
#
#    print 'Computed vanishing point: '
#    print np.concatenate( (vps,tp([conf])) , 1 )
#    print 'Total length: %.1f' % (np.sum(lengths))

#    testvp = geom.normalrows(tvp[0,:])
#    print np.inner(lineqs[-10:],testvp)