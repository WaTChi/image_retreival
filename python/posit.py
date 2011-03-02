#!/usr/bin/env python

import cv
import earthMine as em
import pixels
import tags as tg
import render_tags
import cloud


db = tg.TagCollection('/media/DATAPART2/Research/app/code/tags.csv')

img = "37.8695551919,-122.266734533-0004sift.txt"
jpg = '/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829/37.8695551919,-122.266734533-0004.jpg'
#img = "37.87274692,-122.268484938-0009sift.txt"
#jpg = '/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829/37.87274692,-122.268484938-0009.jpg'
#img = "37.8695529328,-122.268044238-0010sift.txt"
#jpg = '/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829/37.8695529328,-122.268044238-0010.jpg'
#img = "37.8696156756,-122.266254025-0009sift.txt"
#jpg = '/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829/37.8696156756,-122.266254025-0009.jpg'
#img = "37.8696422624,-122.267852592-0009sift.txt"
#jpg = '/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829/37.8696422624,-122.267852592-0009.jpg'
source = render_tags.EarthmineImageInfo(jpg, jpg[:-4] + '.info')
timg = render_tags.TaggedImage(jpg, source, db)
v = {'view-location':{'lat':timg.lat, 'lon':timg.lon, 'alt': timg.alt}}
tags= timg.get_frustum()
#for t in tags:
#    print t
timg.draw(timg.map_tags_camera(), '/media/DATAPART2/out.png')
FOCAL_LENGTH=timg.focal_length
FOCAL_LENGTH=2000
print "FOCAL LENGTH: {0}".format(FOCAL_LENGTH)

# read in points
px = pixels.PixelMap('/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829')
rawlocs = px.open(img)
#filter out ones w/o 3d points
locsar = filter(lambda x: x[1], rawlocs.items())
def test():



    # change 3d coordinate systems
    c = map(lambda x: ({'x':x[0][0], 'y':x[0][1]}, x[1]), locsar)
    locpairs = em.ddImageLocstoLPT(v, c)
    #translate to POSIT Specs
    pts2d = map(lambda x: x[0], locpairs)
    translation2d = (384,256)
    pts2d = map(lambda x: (x[0][0]-translation2d[0], x[0][1]- translation2d[1]), locpairs)
    translation3d = locpairs[0][1]
    pts3d = map(lambda x: tuple(x[1] - translation3d), locpairs)
    #pts3d = map(lambda x: tuple(x[1]), locpairs)

    #convert tags to correct coordinate system
    c = map(lambda x: ({'x':0, 'y':0}, {'lat':x.lat, 'lon':x.lon, 'alt':x.alt}), tags)
    taglocpairs = em.ddImageLocstoLPT(v, c)
    tagpts3d = map(lambda x: tuple(x[1] - translation3d), taglocpairs)


    positobj = cv.CreatePOSITObject(pts3d)
    rotMat, transVec = cv.POSIT(positobj, pts2d, FOCAL_LENGTH, (cv.CV_TERMCRIT_EPS, 0, 0.000001))


    print "rotation matrix:\t{0}".format(rotMat)
    print "translation matrix:\t{0}".format(transVec)

    #change rotMat to cvMat
    rotMatCV = cv.CreateMat(3, 3, cv.CV_64F)
    for i in range(0, 3):
        for j in range(0, 3):
            cv.Set2D(rotMatCV,i, j, rotMat[i][j])
    #convert rotMatrix to rotVector
    rotVec = cv.CreateMat(3, 1, cv.CV_64F)
    cv.Rodrigues2(rotMatCV, rotVec)

    #change transVec to cvMat
    transVecCV = cv.CreateMat(1, 3, cv.CV_64F)
    for i in range(0, 3):
        transVecCV[0, i] = transVec[i]

    #camera matrix
    cameratrans = cv.CreateMat(3, 3, cv.CV_64F)
    cv.SetIdentity(cameratrans)
    cameratrans[0,0]=FOCAL_LENGTH
    cameratrans[1,1]=FOCAL_LENGTH

    #distCoefs
    distCoef = cv.CreateMat(4, 1, cv.CV_64F)
    cv.SetZero(distCoef)

    #change 3d coordinate data format
    pts3d_mat = cv.CreateMat(len(pts3d), 1, cv.CV_64FC3)
    for i, m in enumerate(pts3d):
        cv.Set2D(pts3d_mat, i, 0, cv.Scalar(*m))

    #project points
    d2 = cv.CreateMat(pts3d_mat.rows, 1, cv.CV_64FC2)
    cv.ProjectPoints2(pts3d_mat, rotVec, transVecCV, cameratrans, distCoef, d2)

    xerrors = []
    yerrors = []
    for i in range(0, d2.rows):
    #    print "project to:\t {0}".format(d2[i,0])
    #    print "2dloc:\t\t {0}".format(pts2d[i])
        xerror=abs(pts2d[i][0]-(d2[i,0][0]))
        yerror=abs(pts2d[i][1]-(d2[i,0][1]))
    #    print "error:\t\t {0}, {1}".format(xerror, yerror)
        xerrors.append(xerror)
        yerrors.append(yerror)
    print "avg xerror:\t {0}".format(sum(xerrors)/len(xerrors))
    print "avg yerror:\t {0}".format(sum(yerrors)/len(yerrors))

    #tag change 3d coordinate data format
    tagpts3d_mat = cv.CreateMat(len(tagpts3d), 1, cv.CV_64FC3)
    for i, m in enumerate(tagpts3d):
        cv.Set2D(tagpts3d_mat, i, 0, cv.Scalar(*m))
    #project points
    d2 = cv.CreateMat(tagpts3d_mat.rows, 1, cv.CV_64FC2)
    cv.ProjectPoints2(tagpts3d_mat, rotVec, transVecCV, cameratrans, distCoef, d2)
    ntags=[]
    for i in range(0, d2.rows):
    #    print "project to:\t {0}".format(d2[i,0])
    #    ntags.append((tags[i], (0, d2[i,0])))
        ntags.append((tags[i], (0, (d2[i,0][0]+translation2d[0],d2[i,0][1]+translation2d[1]))))
    return ntags

#timg.draw(test(), '/media/DATAPART2/out2.png')

cloud.setkey(api_key=2160, api_secretkey='d3497353fc98fc4f3d62561c925c97ecd910cfbb')
jid = cloud.call(test) #a jid identifies your job (a function)
timg.draw(cloud.result(jid), '/media/DATAPART2/out2.png')
