#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="aaronh"
__date__ ="$Sep 30, 2011 2:58:41 PM$"

import extract_oakland
import os

if __name__ == "__main__":

    panopath = '/media/DATAPART2/Research/collected_images/earthmine-oakland-plane/oakland_sphericals/1000006383548'
    outpath = os.path.join(panopath,'sample')
    orientation = [166.172771,1.309877,1.37411]
    view = 3
    detail = 1
    extract_oakland.genView(panopath,outpath,orientation,view,detail)