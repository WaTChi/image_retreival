# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="tyler"
__date__ ="Oct 26, 2011"

import Image
import time
import numpy as np
import numpy.linalg as alg
import numpy.random as rnd
import scipy as sp
import scipy.signal as sig
import lsd
import math
import random

def vanishing_points(imgpath,Kcal='default'):
    
    """
    %% compile
    mex -O -output lsd lsd_matlab.c lsd.c
    %% test
    clear all; close all; clc;
    chairs = imread('SnowCity.jpeg');
    imshow(chairs); hold on;
    tic
    lines = lsd(double(chairs));
    %%lines = LineSegmentDetection(double(chairs), 0.8, 0.6, 2.0, 22.5, 0.0, 0.7, 1024, 255.0, 0.0);
    %%disp(lines);
    t = toc;
    disp(['[lsd] ',num2str(t),' seconds elapsed.']);
    nl = size(lines,1);
    disp(nl);
    for i=1:nl
        plot(lines(i,1:2:4),lines(i,2:2:4),'r-');
    end
    """
    
    # takes a long time

    # parameters
    #Gm_threshold = .2
    #vptol = .25
    #maxvps = 20
    #maxlines = 1000
    
    """
    extern void free_image_double(image_double i);
    extern void free_ntuple_list(ntuple_list in);

    extern image_double makeImageDoubleForPython(unsigned int xsize, unsigned int ysize, double *data);
    extern image_double new_image_double(unsigned int xsize, unsigned int ysize);

    extern ntuple_list lsd(image_double image);
    extern ntuple_list lsd_scale(image_double image, double scale);
    """

    # get BW matrix from imgpath
    timeStart = time.time()
    img = np.asarray( Image.open(imgpath).convert('L') ).astype(np.float64) / 256   # L is for black
    print "conversion time is ", time.time() - timeStart
    print "image data type is ", img.dtype
    
    width = img.shape[0]
    length = img.shape[1]
    
    print "image size is ", img.shape
            
    imageAsDoubleArray = lsd.new_image_double_s(width*length)
        
    start = time.time()   
    
    lsd.TheVector_
            
    #flattens to 1-D array, waste of time?
    for i in xrange(width*length):
        offsetX = i%width
        offsetY = math.floor(i/width)
        imageAsDoubleArray[i] = img[offsetX, offsetY]
        if i%1000000 == 0:
            print "i is ", i, " value is ", img[offsetX, offsetY]
            
    imageDouble = lsd.image_double_s(imageAsDoubleArray, width, length)
    
    print "Flatten time is ", time.time() - start            
                        
    answer = lsd.lsd(imageDouble)
    
    values = lsd.ntuple_list_s_values_get()
    
    answer.
    
    #lsd.free_image_double(imageFromC)
    
    #for value in answer:
     #   print value
                    
    #~doublePointer(imageAsDoubleArray)
    
    """  
    
    # Generate unset inputs
    if type(Kcal) is str:
        Kcal = np.array([[ img.shape[1] , 0 , img.shape[1]/2. ],
                         [ 0 , img.shape[0] , img.shape[0]/2. ],
                         [ 0 ,            0 ,              1. ]])

    # compute inverse camera matrix
    Kinv = alg.inv(Kcal)

    # filters, detecting edges
    Hx = np.array([[-1,-2,-1],[0,0,0],[1,2,1]])/8.0
    Hy = Hx.transpose()

    # x, y line detected images
    Gx = sig.convolve2d(img,Hx,'same')
    Gy = sig.convolve2d(img,Hy,'same')
    Gm = np.sqrt( Gx**2 + Gy**2 )

    # Indices with sufficient gradient
    idx = np.nonzero( Gm > Gm_threshold )
    pix = np.transpose( idx )

    # pull out indices that matter and compute 3d rays of pixels
    Gx = Gx[idx[0],idx[1]]
    Gy = Gy[idx[0],idx[1]]
    Gm = Gm[idx[0],idx[1]]
    npts = len(Gm)
    rays = np.dot( Kinv , np.array( [ pix[:,1] , pix[:,0] , np.ones(npts) ] ) ).transpose()  # ray from camera center to the point

    # choose maxlines random vanishing lines if there are too many, randomly
    if npts > maxlines:
        rnd_idx = rnd.permutation(range(pix.shape[0]))[:maxlines]
        Gx = Gx[rnd_idx]
        Gy = Gy[rnd_idx]
        Gm = Gm[rnd_idx]
        rays = rays[rnd_idx,:]

    # number of vanishing lines
    npts = len(Gm)
    rngN = range(npts)
    
    # compute vanishing line equations and their normalized versions
	# Potential Improvement: by eliminating vertical outliers at this point
    lines = np.array( [ -Gy , Gx , Gy*rays[:,0]-Gx*rays[:,1] ]  ).transpose()
    lines = lines / np.tile( (np.sum(lines**2,1)**0.5)[:,np.newaxis] , [1,3] )

    # RANSAC initialization
    bvps = np.array([[np.inf,np.inf,np.inf] for i in xrange(maxvps)])  	# bestVanishingPoints
    berrs = np.array([np.inf for i in xrange(maxvps)])			# errors
    bnumi = np.zeros(maxvps)						# number of inliers
    iter, maxiter, maxerr, minfrac = 0, 1000, .03, .05
    stoperr, minfit = 0.1*np.sqrt(maxerr)/(minfrac*maxlines), max(3,int(minfrac*npts))

    # RANSAC loop
    while iter < maxiter:
        iter += 1
        i0, i1 = tuple(rnd.permutation(rngN)[:2])
        line0, line1 = lines[i0,:], lines[i1,:]
        mvp = np.cross(line0,line1)  # maybe vanishing point
		# Potential Improvement: if this is facing up or somehow not correct we can skip the rest?
        errs = np.abs( np.dot( lines , mvp ) )
        imask = errs < maxerr
        numi = sum(imask)
        if numi > minfit:
            ilines = lines[imask,:]
            [u,s,vt] = alg.svd( ilines )
            ivp = vt[2,:]
            ierrs = np.abs( np.dot( ilines , ivp ) )
            center_bias = alg.norm(ivp[:2])
            ierr = np.mean(ierrs) / (center_bias*numi)
            if ierr < berrs[-1]:
                # check to see if this is same as another vp
                sinVP = np.sqrt( 1 - np.dot( bvps , ivp )**2 )
                if (sinVP<vptol).any():
                    idx_replace = np.argmin(sinVP)
                    if ierr < berrs[idx_replace]:
                        berrs[idx_replace] = ierr
                        bvps[idx_replace,:] = ivp
                        bnumi[idx_replace] = numi
                        shuffle = np.argsort(berrs)
                        berrs = berrs[shuffle]
                        bvps = bvps[shuffle,:]
                        bnumi = bnumi[shuffle]
                else:
                    berrs[-1] = ierr
                    bvps[-1,:] = ivp
                    bnumi[-1] = numi
                    shuffle = np.argsort(berrs)
                    berrs = berrs[shuffle]
                    bvps = bvps[shuffle,:]
                    bnumi = bnumi[shuffle]
                if (berrs*bnumi < stoperr).all():
                    print 'Stopped early.'
                    break

    # remove indices without vanishing points
    valid_idx = ( bnumi > 0 )
    bvps = bvps[valid_idx,:]
    berrs = berrs[valid_idx]
    bnumi = bnumi[valid_idx]

    # remove overlapping vanishing points
    iter = 0
    while iter < len(bnumi)-1:
        vp = bvps[iter]
        sinVP = 1 - np.sqrt( np.dot( bvps[iter+1:] , ivp )**2 )
        idx_replace = np.nonzero( sinVP < vptol )[0]
        if (sinVP<vptol).any():
            idx_delete = np.argmin(sinVP)+1
            bnumi = np.delete(bnumi,idx_delete)
            berrs = np.delete(berrs,idx_delete)
            bvps = np.delete(bvps,idx_delete,0)
        iter += 1

    print 'Extracting vanishing points took %.1f seconds.' % (time.time()-start)
    
    """
    
    return None, None, None

    #return bvps, berrs, bnumi
