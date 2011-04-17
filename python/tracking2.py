import sys
import cv
import time
import numpy as np
import os.path
import os
import client
import info
import bundlerWrapper

BAD_HOMOGRAPHY_DET_THRESHOLD = .01

# parameters for cv.GoodFeaturesToTrack
MAX_CORNERS = 500
QUALITY = 0.01
MIN_DISTANCE = 5.0
BLOCK_SIZE = 3
HARRIS_PARAM = 0.04

# parameters for cv.FindCornerSubPix
WIN = (20, 20)
REFINE_NUM_ITERS = 20
REFINE_ACCURACY = 0.03
REFINE_CRITERIA = (cv.CV_TERMCRIT_ITER | cv.CV_TERMCRIT_EPS,
                   REFINE_NUM_ITERS,
                   REFINE_ACCURACY)

# parameters for cv.CalcOpticalFlowPyrLK
FLOW_NUM_ITERS = 20
FLOW_ACCURACY = 0.03
FLOW_CRITERIA = (cv.CV_TERMCRIT_ITER | cv.CV_TERMCRIT_EPS,
                 FLOW_NUM_ITERS,
                 FLOW_ACCURACY)
FLAGS = 0

# parameters for cv.FindHomography
METHOD = cv.CV_RANSAC
THRESHOLD = 9
STATUS = None


PYRAMIDS = 7
VERBOSE = False
#FILESET = ("/home/etzeng/tracking/img0.jpg","/home/etzeng/tracking/img1.jpg","/home/etzeng/tracking/test/OpticalFlow1.jpg")
#FILESET = ("/home/etzeng/tracking/img02.jpg","/home/etzeng/tracking/img03.jpg","/home/etzeng/tracking/test/OpticalFlow3.jpg")
#FILESET = ("/home/etzeng/tracking/img02s.jpg","/home/etzeng/tracking/img03s.jpg","/home/etzeng/tracking/img03s.jpg")
#FILESET = ("/home/etzeng/tracking/imga.jpg","/home/etzeng/tracking/imgb.jpg","/home/etzeng/tracking/test/OpticalFlowB.jpg")
#FILESET = ("/home/etzeng/tracking/img04.jpg","/home/etzeng/tracking/img05.jpg","/home/etzeng/tracking/test/img05.jpg")
#FILESET = ("/home/etzeng/tracking2/img00.jpg","/home/etzeng/tracking2/img01.jpg","/home/etzeng/tracking2/img01.jpg") # GOOD
#FILESET = ("/home/etzeng/tracking3/100_0437.jpg", "/home/etzeng/tracking3/100_0438.jpg", "/home/etzeng/tracking3/100_0438.jpg")
#FILESET = ("/home/etzeng/tracking4/100B0450.jpg", "/home/etzeng/tracking4/100B0451.jpg", "/home/etzeng/tracking4/100B0451.jpg")
DEFAULT_FILESET = ("/home/etzeng/tracking5/IMG_4571.jpg", "/home/etzeng/tracking5/IMG_4572.jpg", "/home/etzeng/tracking5/IMG_4572.jpg")

   
def find_features(img, num_corners=MAX_CORNERS, quality=QUALITY,
                  min_dist=MIN_DISTANCE, mask=None, block_size=BLOCK_SIZE,
                  use_harris=False, harris_param=HARRIS_PARAM):
    image_size = cv.GetSize(img)
    # cv.GoodFeaturesToTrack needs two temp images
    eig_image = cv.CreateImage(image_size, cv.IPL_DEPTH_32F, 1)
    tmp_image = cv.CreateImage(image_size, cv.IPL_DEPTH_32F, 1)
    return cv.GoodFeaturesToTrack(img, eig_image, tmp_image, num_corners,
                                  quality, min_dist, mask, block_size, 
                                  1 if use_harris else 0, harris_param)

def refine_features(img, corners, win=WIN, zero_zone=(-1, -1),
                    criteria=REFINE_CRITERIA):
    return cv.FindCornerSubPix(img, corners, win, zero_zone, criteria)

def compare_pair(file1, file2):
    start = time.time()
    img1 = cv.LoadImage(file1, cv.CV_LOAD_IMAGE_GRAYSCALE)
    img2 = cv.LoadImage(file2, cv.CV_LOAD_IMAGE_GRAYSCALE)
    image_size = cv.GetSize(img1)
    # need another image to draw the flow vectors onto, so grab the second image again
    img_flow = cv.LoadImage(file2, cv.CV_LOAD_IMAGE_UNCHANGED)
    corners1 = find_features(img1)
    corners1 = refine_features(img1, corners1)
    # this step is the magical one
    corners2, features_found, feature_errors = lucas_kanade(img1, img2, corners1)
    homography = find_homography(corners1, corners2)
    # bundler stuff
    r, t = bundlerWrapper.find_r_and_t(file1, file2)
    if r == None and t == None:
        print "bad pair"
        return None, None
    print "==="
    print "rotation:", r
    print "translation:", t
    # homographyDecomposition.extract_R_and_T(homography)
    # decomposition = decompose_projection(homography)
    # print decomposition
    print "TIME ELAPSED: ", str(time.time() - start)
    img1 = cv.LoadImage(file1, cv.CV_LOAD_IMAGE_UNCHANGED)
    img2 = cv.LoadImage(file2, cv.CV_LOAD_IMAGE_UNCHANGED)
    display_tracking(image_size, img1, img2, img_flow, homography)
    return r, t
    
def lucas_kanade(img1, img2, corners1, win=WIN,
                 num_pyramids=PYRAMIDS, criteria=FLOW_CRITERIA, flags=FLAGS):
    pyr_size = (img1.width + 8, img2.height / 3) # magic formula provided by openCV
    pyr1 = cv.CreateImage(pyr_size, cv.IPL_DEPTH_32F, 1)
    pyr2 = cv.CreateImage(pyr_size, cv.IPL_DEPTH_32F, 1)
    return cv.CalcOpticalFlowPyrLK(img1, img2, pyr1, pyr2, corners1, win,
                                   num_pyramids, criteria, flags)

def find_homography(corners1, corners2, method=METHOD, threshold=THRESHOLD,
                    status=STATUS):
    array1 = np.mat(corners1)
    array2 = np.mat(corners2)
    homography = cv.CreateMat(3, 3, cv.CV_64F)
    cv.FindHomography(cv.fromarray(array1), cv.fromarray(array2), homography,
                      METHOD, THRESHOLD)
    if VERBOSE:
        # TODO: find a less obnoxious way to do this
        print "Homography matrix:"
        for i in range(3):
            print homography[i,0],
            print homography[i,1],
            print homography[i,2]
    #find_svd(homography)
    return homography

def find_svd(matrix):
    W = cv.CreateMat(3, 3, cv.CV_64F)
    U = cv.CreateMat(3, 3, cv.CV_64F)
    V = cv.CreateMat(3, 3, cv.CV_64F)
    cv.SVD(matrix, W, U, V, cv.CV_SVD_U_T | cv.CV_SVD_V_T) 
    print "W:"
    for i in range(3):
        print W[i,0], W[i,1], W[i,2]
    print "U:"
    for i in range(3):
        print U[i,0], U[i,1], U[i,2]
    print "V:"
    for i in range(3):
        print V[i,0], V[i,1], V[i,2]
    return W, U, V

def decompose_projection(matrix):
    H = cv.CreateMat(3, 4, cv.CV_64F)
    for i in range(3):
        for j in range(3):
            H[i,j] = matrix[i,j]
        H[i,3] = 0
    K = cv.CreateMat(3, 3, cv.CV_64F)
    R = cv.CreateMat(3, 3, cv.CV_64F)
    T = cv.CreateMat(4, 1, cv.CV_64F)
    euler = cv.DecomposeProjectionMatrix(H, K, R, T)
    print "trans:", K[0,0], K[0,1], K[0,2]
    return K, R, T, euler
 
def display_tracking(img_size, img1, img2, img_flow, homography):
    for x in range(0, img_size[0], 100):
        for y in range(0, img_size[1], 100):
            cv.Circle (img1, (x, y), 3, (0, 255, 0, 0), -1, 8, 0)
            point = cv.CreateMat(3, 1, cv.CV_64F)
            point[0,0] = x
            point[1,0] = y
            point[2,0] = 1
            newpoint = cv.CreateMat(3, 1, cv.CV_64F)
            cv.MatMul(homography, point, newpoint)
            cv.Circle (img2, (int(newpoint[0,0]), int(newpoint[1,0])), 3, (0, 255, 0, 0), -1, 8, 0)
            cv.Line(img_flow, (x, y), (int(newpoint[0,0]),int(newpoint[1,0])), cv.CV_RGB(255,0,0), 2)
    cv.Rectangle(img_flow, (0,0), (150,25), cv.CV_RGB(0,0,0), thickness=-1)
    cv.PutText(img_flow, "Good?: " + str(isHomographyGood(homography)), (0, 20), cv.InitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.75, 0.75, thickness=2), cv.CV_RGB(255,255,255))
    cv.NamedWindow("Image1", 0)
    cv.NamedWindow("Image2", 0)
    cv.NamedWindow("Flow",0)
    cv.ResizeWindow("Image1",int(0.5*img_size[0]), int(0.5*img_size[1]))
    cv.ResizeWindow("Image2",int(0.5*img_size[0]), int(0.5*img_size[1]))
    cv.ResizeWindow("Flow",int(0.5*img_size[0]), int(0.5*img_size[1]))
    cv.ShowImage("Image1", img1)
    cv.ShowImage("Image2", img2)
    cv.ShowImage("Flow", img_flow)
    cv.WaitKey(0)

def compare_dir(dirname):
    #print os.listdir(dirname)
    #lat, lon = read_coords(dirname)
    files = sorted([os.path.join(dirname, file) for file in os.listdir(dirname)
             if os.path.splitext(file)[1] in (".jpg",".JPG")])
    #matchedimg, matches = client.match(os.path.splitext(files[0])[0] + "sift.txt", lat, lon)
    #print matchedimg, matches
    #print files
    if len(files) < 2:
        print "Not enough images in directory!"
    else:
        total_r = [0, 0, 0]
        total_t = [0, 0, 0]
        while True:
            r, t = compare_pair(files[0], files[1])
            if r == None and t == None:
                files[1:] = files[2:]
                continue
            for i in range(3):
                total_r[i] += r[i]
                total_t[i] += t[i]
            print "cumulative r:", total_r
            print "cumulative t:", total_t
            files = files[1:]
    
def read_coords(dirname):
    # TODO: Pull coordinates out of filename if coords.txt doesn't exist
    coords_file = os.path.join(dirname, "coords.txt")
    if os.path.exists(coords_file):
        lines = open(coords_file).readlines()
        return float(lines[1]), float(lines[2])


def disparity(file1, file2):
    sc = cv.CreateStereoBMState(preset=cv.CV_STEREO_BM_BASIC, numberOfDisparities=0)
    #sc2 = cv.CreateStereoGCState(16, 2)
    img1 = cv.LoadImage(file1, cv.CV_LOAD_IMAGE_GRAYSCALE)
    img2 = cv.LoadImage(file2, cv.CV_LOAD_IMAGE_GRAYSCALE)
    img_sz = cv.GetSize(img1)
    img_out = cv.CreateImage(img_sz, cv.IPL_DEPTH_32F, 1)
    img_out2 = cv.CreateImage(img_sz, cv.IPL_DEPTH_32F, 1)
    img_out3 = cv.CreateImage(img_sz, cv.IPL_DEPTH_32F, 1)
    cv.FindStereoCorrespondenceBM(img1, img2, img_out, sc)
    #cv.FindStereoCorrespondenceGC(img1, img2, img_out2, img_out3, sc2)
    cv.NamedWindow("Image1", 0)
    cv.NamedWindow("Image2", 0)
    cv.NamedWindow("Disparity",0)
    cv.ResizeWindow("Image1",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
    cv.ResizeWindow("Image2",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
    cv.ResizeWindow("Disparity",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
    cv.ShowImage("Image1", img1)
    cv.ShowImage("Image2", img2)
    cv.ShowImage("Disparity", img_out)
    cv.WaitKey(0)
 #    cv.NamedWindow("Image1", 0)
 #   cv.NamedWindow("Image2", 0)
 #   cv.NamedWindow("Disparity1",0)
 #   cv.NamedWindow("Disparity2",0)
 #   cv.ResizeWindow("Image1",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
 #   cv.ResizeWindow("Image2",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
 #   cv.ResizeWindow("Disparity1",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
 #   cv.ResizeWindow("Disparity2",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
 #   cv.ShowImage("Image1", img1)
 #   cv.ShowImage("Image2", img2)
 #   cv.ShowImage("Disparity1", img_out2)
 #   cv.ShowImage("Disparity2", img_out3)
 #   cv.WaitKey(0)
 

def disparity_dir(dirname):
    files = sorted([os.path.join(dirname, file) for file in os.listdir(dirname)
                    if os.path.splitext(file)[1] in (".jpg",".JPG")])
    if len(files) < 2:
        print "Not enough images in directory!"
    else:
        while True:
            disparity(files[0], files[1])
            files.append(files[0])
            files = files[1:]

def isHomographyGood(H):
  pts = [(384,256), (384,361), (489,256)]
  det = np.linalg.det(H)
  if det < BAD_HOMOGRAPHY_DET_THRESHOLD:
    return False
  scale = 1
  dests = []
  for (y,x) in pts:
    dest = H*np.matrix([x,y,1]).transpose()
    try:
      dest = tuple(map(int, (dest[0].item()/dest[2].item(), dest[1].item()/dest[2].item())))
    except ZeroDivisionError:
      dest = (0,0)
    except ValueError:
      dest = (0,0)
    dest = (dest[1]*scale, dest[0]*scale)
    dests.append(dest)
  v1len = np.sqrt((dests[0][0] - dests[1][0])**2 + (dests[0][1] - dests[1][1])**2)
  v2len = np.sqrt((dests[0][0] - dests[2][0])**2 + (dests[0][1] - dests[2][1])**2)
  length_ratio = v1len / v2len
  if length_ratio > 5 or 1/length_ratio > 5:
    return False
  a = (dests[1][0] - dests[0][0], dests[1][1] - dests[0][1])
  b = (dests[2][0] - dests[0][0], dests[2][1] - dests[0][1])
  aLen = (a[0]**2 + a[1]**2)**.5
  bLen = (b[0]**2 + b[1]**2)**.5
  angle_delta = np.arccos(np.dot(a, b)/(aLen*bLen))*180.0/np.pi
  if angle_delta < 30 or angle_delta > 160:
    return False
  return True

if __name__ == '__main__':
    if len(sys.argv) == 1:
        # use DEFAULT_FILESET
        compare_pair(DEFAULT_FILESET[0], DEFAULT_FILESET[1])
        disparity(DEFAULT_FILESET[0], DEFAULT_FILESET[1])
    elif len(sys.argv) == 2 and os.path.isdir(sys.argv[1]):
        # go pairwise through directory
        compare_dir(sys.argv[1])
    elif len(sys.argv) == 3:
        # compare two images
        disparity_dir(sys.argv[1])
 
