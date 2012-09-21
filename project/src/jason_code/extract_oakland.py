import Image
import numpy as np
import scipy.interpolate as spint
import numpy.linalg as alg
import geom
import os
import time
import cv2.cv as cv

if __name__ == "__main__":

    panopath = '/media/DATAPART2/Research/collected_images/earthmine-oakland-plane/oakland_sphericals/1000006374536'
    outpath = panopath + '/test'
    orientation = [108.713523,1.544953,0.888563]
    view = 10
    detail = 1
    #genView(panopath,outpath,orientation,view,detail)
    rasterpath = panopath + '/raster.jpg'
#    depthpath = panopath + '/depth.png'
    planepath = panopath + '/plane_pano.png'
    raster = Image.open(rasterpath)
#    depth = Image.open(depthpath)
    plane = Image.open(planepath)
#    test = Image.new('RGB',(2500,1200))
#    test.save(outpath + '/test.jpg','jpeg')
    pix = np.array(raster.load())
    print pix
#    print raster.mode
    print pix[:,498:502,248:252]

def plane_idx(panopath):
    return


def genView(panopath, outpath, orientation, viewdir, detail):
    # panopath: location of raster, depth and plane_pano images
    # outpath: location to put generated view, depth, and plane images
    # orientation: [yaw, pitch, roll] of panorama (in degrees)
    # viewdir: clock-based view; in set [2,3,4,8,9,10] for database
    # detail: 0 = 768x512 with 60d fov, 1 = 2500x1200 with 90d fov

    # local constants
    start = time.time()
    pi = np.pi

    width, height, fov = (2500, 1200, 90) if detail else (768, 512, 60)

    # view details
    Rpano = geom.RfromYPR(orientation[0],orientation[1],orientation[2])
    Yview = np.mod( orientation[0] + 30*viewdir, 360 )
    Rview = geom.RfromYPR(Yview, 0, 0)
    Kview = geom.cameramat(width, height, fov*pi/180)
    Kinv = alg.inv(Kview)

    # Load image pano, depth pano, and plane pano images
    cvIP = cv.LoadImageM( os.path.join(panopath,'raster.jpg'), cv.CV_LOAD_IMAGE_UNCHANGED )
    cvDP = cv.fromarray( np.asarray( Image.open( os.path.join(panopath,'depth.png') ) ) )
    pp = np.asarray( Image.open( os.path.join(panopath,'plane_pano.png') ) ).copy()
    vals = set(list(pp.reshape(-1)))
    vals.remove(255)
    gnd = max(vals)
    pp[pp==gnd] = np.uint8(0)
    cvPP =  cv.fromarray(pp)

    # load pixel map
    pix = np.append(np.array(np.meshgrid(range(width),range(height))).reshape(2,-1),np.ones([1,width*height]),0)

    midpoint = time.time()
    print 'Loading pano images took ' + str(midpoint-start) + ' seconds.'

    # Create output openCV matrices
    cvI = cv.CreateMat(height,width,cv.CV_8UC3)
    cvD = cv.CreateMat(height,width,cv.CV_32SC1)
    cvP = cv.CreateMat(height,width,cv.CV_8UC1)

    # compute mappings
    ray = np.dot( np.transpose(Rpano), np.dot( Rview, np.dot( Kinv, pix ) ) )
    yaw, pitch = np.arctan2( ray[0,:] , ray[2,:] ) , np.arctan2( -ray[1,:] , np.sqrt((np.array([ray[0,:],ray[2,:]])**2).sum(0)) )
    ix, iy = cv.fromarray(np.array(8192/(2*pi)*(pi+yaw),np.float32).reshape(height,width)), cv.fromarray(np.array(4096/pi*(pi/2-pitch),np.float32).reshape(height,width))
    dx, dy = cv.fromarray(np.array(5000/(2*pi)*(pi+yaw),np.float32).reshape(height,width)), cv.fromarray(np.array(2500/pi*(pi/2-pitch),np.float32).reshape(height,width))
    px, py = cv.fromarray(np.array(1000/(2*pi)*(pi+yaw),np.float32).reshape(height,width)), cv.fromarray(np.array( 500/pi*(pi/2-pitch),np.float32).reshape(height,width))

    # call remap function
    cv.Remap(cvIP,cvI,ix,iy,cv.CV_INTER_CUBIC)
    cv.Remap(cvDP,cvD,dx,dy,cv.CV_INTER_NN)
    cv.Remap(cvPP,cvP,px,py,cv.CV_INTER_NN)

    # write images to file
    Image.fromarray(np.array(cvI)[:,:,[2,1,0]]).save(os.path.join(outpath,'view.jpg'),'jpeg')
    Image.fromarray(np.array(cvD)).save(os.path.join(outpath,'depth.png'),'png')
    Image.fromarray(np.array(cvP)).save(os.path.join(outpath,'plane.png'),'png')

    print 'Generating views from pano took ' + str(time.time()-midpoint) + ' seconds.'