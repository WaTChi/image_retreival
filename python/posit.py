import os
import os.path
#!/usr/bin/env python

import cv
import earthMine as em
import pixels
import tags as tg
import render_tags
import cloud
import util



def test(locsar, tags, viewpoint):

    # change 3d coordinate systems
    c = map(lambda entry: ({'x':entry[0][0], 'y':entry[0][1]}, entry[1]), locsar)
    locpairs = em.ddImageLocstoLPT(viewpoint, c)
    #translate to POSIT Specs
    pts2d = map(lambda x: x[0], locpairs)
    #TODO: constant here
    translation2d = (384,256)
    pts2d = map(lambda x: (x[0][0]-translation2d[0], x[0][1]- translation2d[1]), locpairs)
    translation3d = locpairs[0][1]
    pts3d = map(lambda x: tuple(x[1] - translation3d), locpairs)

    #convert tag coordinate system
    c = map(lambda x: ({'x':0, 'y':0}, {'lat':x.lat, 'lon':x.lon, 'alt':x.alt}), tags)
    taglocpairs = em.ddImageLocstoLPT(viewpoint, c)
    tagpts3d = map(lambda x: tuple(x[1] - translation3d), taglocpairs)

    #compute POSIT
    positobj = cv.CreatePOSITObject(pts3d)
    rotMat, transVec = cv.POSIT(positobj, pts2d, FOCAL_LENGTH, (cv.CV_TERMCRIT_EPS, 0, 0.000001))

#    print "rotation matrix:\t{0}".format(rotMat)
#    print "translation matrix:\t{0}".format(transVec)

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

    #compute self errors
    xerrors = []
    yerrors = []
    for i in range(0, d2.rows):
        xerror=abs(pts2d[i][0]-(d2[i,0][0]))
        yerror=abs(pts2d[i][1]-(d2[i,0][1]))
        xerrors.append(xerror)
        yerrors.append(yerror)
    print "avg xerror:\t {0}".format(sum(xerrors)/len(xerrors))
    print "avg yerror:\t {0}".format(sum(yerrors)/len(yerrors))


    #change tag 3d coordinate data format
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

#cloud.setkey(api_key=2160, api_secretkey='d3497353fc98fc4f3d62561c925c97ecd910cfbb')
#jid = cloud.call(test) #a jid identifies your job (a function)
#timg.draw(cloud.result(jid), '/media/DATAPART2/out2.png')


db = tg.TagCollection('/media/DATAPART2/Research/app/code/tags.csv')
imgdir='/media/DATAPART2/Research/cells/g=100,r=d=236.6/37.8714488812,-122.266998471'
infodir='/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829'
for file in util.getJPGFileNames(imgdir)[0:100]:
    jpg = os.path.join(imgdir, file)
    img=file[:-4]+'sift.txt'
    info=os.path.join(infodir, file[:-4]+'.info')
    out=file[:-4]+'tagged.png'
    out2=file[:-4]+'gttagged.png'

    if not (os.path.exists(info) and os.path.exists(jpg)):
        continue
    try:
        source = render_tags.EarthmineImageInfo(jpg, info)
        timg = render_tags.TaggedImage(jpg, source, db)
    except:
        continue
    v = {'view-location':{'lat':timg.lat, 'lon':timg.lon, 'alt': timg.alt}}
    tags= timg.get_frustum()
    
    try:
        timg.draw(timg.map_tags_camera(), os.path.join('/media/DATAPART2/jz/posit3/',out2))
    except:
        print "error"
    FOCAL_LENGTH=timg.focal_length
    #FOCAL_LENGTH=2000
    print "FOCAL LENGTH: {0}".format(FOCAL_LENGTH)

    # read in points
    px = pixels.PixelMap('/media/DATAPART2/Research/collected_images/earthmine-fa10.1/37.871955,-122.270829')
    rawlocs = px.open(img)
    #filter out ones w/o 3d points
    locsar = filter(lambda x: x[1], rawlocs.items())

	print locsar[0]
	assert false
	
    print "num 3d pts: {0}".format(len(locsar))
    print "num tags in frustrum: {0}".format(len(tags))

    if len(locsar)>0 and len(tags)>5:
        try:
           timg.draw(test(locsar, tags, v), os.path.join('/media/DATAPART2/jz/posit3/',out))
        except:
            print "error"