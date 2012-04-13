#!/usr/bin/env python

import os
import numpy as np
from multiprocessing import Pool, cpu_count
import numpy.linalg as alg
import cv
import Image
import geom

def extract_panorama(panopath, outbase, panoinfo, detail):
    """Generates raster, plane, and depth views at rot_degrees"""

    print "Processing panorama " + '%d' % panoinfo['pano'][0]

    # panopath: location of raster, depth and plane_pano images
    # outbase: base name to put generated view, depth, and plane images
    # panoinfo: Contains information about the panoramic scene
    # detail: 0 = 768x512 with 60d fov, 1 = 2500x1200 with 90d fov

    # local constants
    pi = np.pi
    width, height, fov = (2500, 1200, 90) if detail else (768, 512, 60)
    
    # pano and view details details
    orientation = [panoinfo['yaw'][0],panoinfo['pitch'][0],panoinfo['roll'][0]]
    Rpano = geom.RfromYPR(orientation[0],orientation[1],orientation[2])
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

    # Create output openCV matrices
    cvI = cv.CreateMat(height,width,cv.CV_8UC3)
    cvD = cv.CreateMat(height,width,cv.CV_32SC1)
    cvP = cv.CreateMat(height,width,cv.CV_8UC1)

    for viewdir in [2,3,4,8,9,10]:
        
        # add to base name and generate info file
        viewname = outbase + '%04d' % viewdir
        gen_info_file(panoinfo, viewname + '.info', detail, 30*viewdir)

        # generate view orientation
        Yview = np.mod( orientation[0] + 30*viewdir, 360 )
        Rview = geom.RfromYPR(Yview, 0, 0)

        # compute mappings
        ray = np.dot( np.transpose(Rpano), np.dot( Rview, np.dot( Kinv, pix ) ) )
        yaw, pitch = np.arctan2( ray[0,:] , ray[2,:] ) , np.arctan2( -ray[1,:] , np.sqrt((np.array([ray[0,:],ray[2,:]])**2).sum(0)) )
        ix, iy = cv.fromarray(np.array(8192/(2*pi)*(pi+yaw),np.float32).reshape(height,width)), cv.fromarray(np.array(4096/pi*(pi/2-pitch),np.float32).reshape(height,width))
        dx, dy = cv.fromarray(np.array(5000/(2*pi)*(pi+yaw),np.float32).reshape(height,width)), cv.fromarray(np.array(2500/pi*(pi/2-pitch),np.float32).reshape(height,width))
        px, py = cv.fromarray(np.array(1000/(2*pi)*(pi+yaw),np.float32).reshape(height,width)), cv.fromarray(np.array( 500/pi*(pi/2-pitch),np.float32).reshape(height,width))

        # call remap function
        cv.Remap(cvIP,cvI,ix,iy,cv.CV_INTER_CUBIC) # if detail else cv.Remap(cvIP,cvI,ix,iy,cv.CV_INTER_AREA)
        cv.Remap(cvDP,cvD,dx,dy,cv.CV_INTER_NN)
        cv.Remap(cvPP,cvP,px,py,cv.CV_INTER_NN)

        # write images to file
        Image.fromarray(np.array(cvI)[:,:,[2,1,0]]).save(viewname+'.jpg','jpeg')
        Image.fromarray(np.array(cvD)).save(viewname+'-depth.png','png')
        Image.fromarray(np.array(cvP)).save(viewname+'-planes.png','png')

def load_meta(spreadsheet):
    """Returns map<panorama_id => record-array"""

    return np.genfromtxt(spreadsheet, dtype=[
        ('pano', '<i8'),
        ('lon', '<f8'),
        ('lat', '<f8'),
        ('alt', '<f8'),
        ('yaw', '<f8'),
        ('pitch', '<f8'),
        ('roll', '<f8'),
        ('gmt', '|S22'),
        ('blue_ratio', '<f8')
    ], delimiter='|', names=None)

def gen_info_file(data, out, detail, rot_degrees):
    """Writes data + rot_degrees to out, formatted in .info format"""

    if os.path.exists(out):
        return

    dict = {
        'is-known-occluded': False,
        'url': {'href': ''},
        'field-of-view': 30.0,
        'image-size': {
            'width': 300,
            'height': 225,
        },
        'view-direction': {
            'yaw': 0.0,
            'pitch': 0.0,
        },
        'view-location': {
            'lat': 37.875507,
            'lon': -122.264883,
            'alt': 0,
        },
        'id': 'tutorial_example_image',
    }
    with open(out, 'w') as f:
        f.write(str(dict))

def convert_panoramas(metadata, dir, outdir, detail):
    """Given a panorama dir, extract rectangular data"""

    pool = Pool(cpu_count())

    data = metadata[metadata['pano'] == int(os.path.basename(dir))]
    assert len(data) == 1

    lat, lon = data['lat'][0], data['lon'][0]
    outbase = os.path.join(outdir, '%f,%f-' % (lat, lon))
    pool.apply_async(extract_panorama,
        [dir, outbase, data, detail])

    pool.close()
    pool.join()

if __name__ == '__main__':
    PANORAMA_DIR = '/media/DATAPART1/oakland/earthmine/oakland_sphericals'
    OUTPUT_DIR = '/media/DATAPART1/oakland/earthmine/rect'
    SPREADSHEET = '/media/DATAPART1/oakland/earthmine/oakland_spherical_metadata.csv'
    HIRES_VIEWS = True
    if HIRES_VIEWS:
        OUTPUT_DIR = OUTPUT_DIR + '_hires'
    metadata = load_meta(SPREADSHEET)
    for dir in os.listdir(PANORAMA_DIR):
        dir = os.path.join(PANORAMA_DIR, dir)
        convert_panoramas(metadata, dir, OUTPUT_DIR, HIRES_VIEWS)
