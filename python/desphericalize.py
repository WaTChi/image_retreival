#!/usr/bin/env python

import os
import numpy as np
from multiprocessing import Pool, cpu_count

def extract_panorama(raster, depth, planes, rot_radians, out):
    """Generates raster, plane, and depth views at rot_radians"""

    print "desphericalizing rot =", rot_radians, out

    raster_out = out + '.jpg'
    depth_out = out + '-depth.png'
    plane_out = out + '-planes.png'

    done = all(map(os.path.exists, [raster_out, depth_out, plane_out]))
    if not done:
        print "TODO hook into aaron's code"

def load_meta(spreadsheet):
    """Returns map<panorama_id => record-array"""

    return np.genfromtxt(spreadsheet, dtype=[
        ('pano', '<i8'),
        ('lat', '<f8'),
        ('lon', '<f8'),
        ('alt', '<f8'),
        ('yaw', '<f8'),
        ('pitch', '<f8'),
        ('roll', '<f8'),
        ('gmt', '|S22'),
        ('blue_ratio', '<f8')
    ], delimiter='|', names=None)

def gen_info_file(data, out, rot_radians):
    """Writes data + rot_radians to out, formatted in .info format"""

    if os.path.exists(out):
        return

    dict = {
        'is-known-occluded': False,
        'url': {'href': ''},
        'field-of-view': 60.0,
        'image-size': {
            'width': 768,
            'height': 512,
        },
        'view-direction': {
            'yaw': (data['yaw'][0] + rot_radians) % (np.pi*2),
            'pitch': 0.0,
        },
        'view-location': {
            'lat': data['lat'][0],
            'lon': data['lon'][0],
            'alt': data['alt'][0],
        },
        'id': '',
    }
    with open(out, 'w') as f:
        f.write(str(dict))

def convert_panoramas(metadata, dir, outdir):
    """Given a panorama dir, extract rectangular data"""

    pool = Pool(cpu_count())

    data = metadata[metadata['pano'] == int(os.path.basename(dir))]
    assert len(data) == 1

    raster = os.path.join(dir, 'raster.jpg')
    depth = os.path.join(dir, 'depth.png')
    planes = os.path.join(dir, 'plane_pano.png')
    lat, lon = data['lat'][0], data['lon'][0]

    for rot in [2,3,4,8,9,10]:
        out = os.path.join(outdir, '%f,%f-%04d' % (lat, lon, rot))
        rot_radians = rot*np.pi/6
        gen_info_file(data, out + '.info', rot_radians)
        pool.apply_async(extract_panorama,
            [raster, depth, planes, rot_radians, out])

    pool.close()
    pool.join()

if __name__ == '__main__':
    PANORAMA_DIR = '/media/DATAPART2/Research/collected_images/earthmine-oakland-plane/oakland_sphericals'
    OUTPUT_DIR = '/media/DATAPART2/Research/collected_images/earthmine-oakland-plane/oakland_rect'
    SPREADSHEET = '/media/DATAPART2/Research/collected_images/earthmine-oakland-plane/oakland_spherical_metadata.csv'
    metadata = load_meta(SPREADSHEET)
    for dir in os.listdir(PANORAMA_DIR):
        dir = os.path.join(PANORAMA_DIR, dir)
        convert_panoramas(metadata, dir, OUTPUT_DIR)
