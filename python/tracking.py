import cv
import time
import numpy as np

MAX_CORNERS = 500
PYRAMIDS = 7
FILESET = ("/home/etzeng/tracking/img0.jpg","/home/etzeng/tracking/img1.jpg","/home/etzeng/tracking/test/OpticalFlow1.jpg")
#FILESET = ("/home/etzeng/tracking/img02.jpg","/home/etzeng/tracking/img03.jpg","/home/etzeng/tracking/test/OpticalFlow3.jpg")
#FILESET = ("/home/etzeng/tracking/img02s.jpg","/home/etzeng/tracking/img03s.jpg","/home/etzeng/tracking/img03s.jpg")
#FILESET = ("/home/etzeng/tracking/imga.jpg","/home/etzeng/tracking/imgb.jpg","/home/etzeng/tracking/test/OpticalFlowB.jpg")
#FILESET = ("/home/etzeng/tracking/img04.jpg","/home/etzeng/tracking/img05.jpg","/home/etzeng/tracking/test/img05.jpg")
#FILESET = ("/home/etzeng/tracking2/img00.jpg","/home/etzeng/tracking2/img01.jpg","/home/etzeng/tracking2/img01.jpg") # GOOD
#FILESET = ("/home/etzeng/tracking3/100_0437.jpg", "/home/etzeng/tracking3/100_0438.jpg", "/home/etzeng/tracking3/100_0438.jpg")
FILESET = ("/home/etzeng/tracking4/100B0450.jpg", "/home/etzeng/tracking4/100B0451.jpg", "/home/etzeng/tracking4/100B0451.jpg")


if __name__ == '__main__':
    start = time.time()
    imgA = cv.LoadImage(FILESET[0],cv.CV_LOAD_IMAGE_GRAYSCALE)
    imgB = cv.LoadImage(FILESET[1],cv.CV_LOAD_IMAGE_GRAYSCALE)
    
    img_sz = cv.GetSize(imgA)
    win_size = 20
    
    imgC = cv.LoadImage(FILESET[2],cv.CV_LOAD_IMAGE_UNCHANGED)
    
    eig_image = cv.CreateImage(img_sz, cv.IPL_DEPTH_32F, 1)
    tmp_image = cv.CreateImage(img_sz, cv.IPL_DEPTH_32F, 1)
    
    corner_count = MAX_CORNERS
  
    cornersA = cv.GoodFeaturesToTrack(imgA, eig_image, tmp_image, corner_count, 0.01, 5.0, None, 3, 0, 0.04)
    cornersA = cv.FindCornerSubPix(imgA, cornersA, (win_size, win_size), (-1, -1), (cv.CV_TERMCRIT_ITER | cv.CV_TERMCRIT_EPS, 20, 0.03))
  
    # Call the Lucas Kanade algorithm
    features_found = []
    feature_errors = []

    pyr_sz = (imgA.width + 8, imgB.height / 3)
    
    pyrA = cv.CreateImage(pyr_sz, cv.IPL_DEPTH_32F, 1)
    pyrB = cv.CreateImage(pyr_sz, cv.IPL_DEPTH_32F, 1)
    
    (cornersB, features_found, feature_errors) = cv.CalcOpticalFlowPyrLK(imgA, imgB, pyrA, pyrB, cornersA, (win_size, win_size), PYRAMIDS, (cv.CV_TERMCRIT_ITER | cv.CV_TERMCRIT_EPS, 20, .3), 0)
#    vel_sz = ((img_sz[0]-10)/1, (img_sz[1]-10)/1)
#    velX, velY = cv.CreateImage(vel_sz, cv.IPL_DEPTH_32F,1), cv.CreateImage(vel_sz, cv.IPL_DEPTH_32F,1)
#    cv.CalcOpticalFlowBM(imgA, imgB, (10, 10), (1, 1), (20, 20), 0, velX, velY)
    arrayA = np.mat(cornersA)
#    print arrayA
    arrayB = np.mat(cornersB)
#    print arrayB
    homography_matrix = cv.CreateMat(3, 3, cv.CV_64F)
    cv.FindHomography(cv.fromarray(arrayA), cv.fromarray(arrayB), homography_matrix, cv.CV_RANSAC, 9)
    print "Homography matrix:"
    for i in range(3):
        print homography_matrix[i,0], homography_matrix[i,1], homography_matrix[i,2]
    
    W = cv.CreateMat(3, 3, cv.CV_64F)
    U = cv.CreateMat(3, 3, cv.CV_64F)
    V = cv.CreateMat(3, 3, cv.CV_64F)
    cv.SVD(homography_matrix, W, U, V, cv.CV_SVD_U_T | cv.CV_SVD_V_T) 
    print "W:"
    for i in range(3):
        print W[i,0], W[i,1], W[i,2]
    print "U:"
    for i in range(3):
        print U[i,0], U[i,1], U[i,2]
    print "V:"
    for i in range(3):
        print V[i,0], V[i,1], V[i,2]
                 
  
    print "TIME ELAPSED: ", str(time.time() - start)

    for x in range(0,img_sz[0], 100):
        for y in range(0,img_sz[1],100):
            cv.Circle (imgA, (x, y), 3, (0, 255, 0, 0), -1, 8, 0)
            point = cv.CreateMat(3, 1, cv.CV_64F)
            point[0,0] = x
            point[1,0] = y
            point[2,0] = 1
            newpoint = cv.CreateMat(3, 1, cv.CV_64F)
            cv.MatMul(homography_matrix, point, newpoint)
            cv.Circle (imgB, (int(newpoint[0,0]), int(newpoint[1,0])), 3, (0, 255, 0, 0), -1, 8, 0)
            cv.Line(imgC, (x, y), (int(newpoint[0,0]),int(newpoint[1,0])), cv.CV_RGB(255,0,0), 2)

#    for i in range(corner_count):
#        if features_found[i] == 0 or feature_errors[i] > 550:
#            #print "Error is {0}".format(feature_errors[i])
#            continue
#        p0 = (cv.Round(cornersA[i][0]), cv.Round(cornersA[i][1]))
#        p1 = (cv.Round(cornersB[i][0]), cv.Round(cornersB[i][1]))
#        cv.Line(imgC, p0, p1, cv.CV_RGB(255,0,0), 2)

    cv.NamedWindow("ImageA", 0)
    cv.NamedWindow("ImageB", 0)
    cv.NamedWindow("LKpyr_OpticalFlow",0)
#    cv.NamedWindow("VelX", 0)
#    cv.NamedWindow("VelY", 0)

    cv.ResizeWindow("ImageA",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
    cv.ResizeWindow("ImageB",int(0.5*img_sz[0]), int(0.5*img_sz[1]))
    cv.ResizeWindow("LKpyr_OpticalFlow",int(0.5*img_sz[0]), int(0.5*img_sz[1]))

    cv.ShowImage("ImageA", imgA)
    cv.ShowImage("ImageB", imgB)
    cv.ShowImage("LKpyr_OpticalFlow", imgC)
#    cv.ShowImage("VelX", velX)
#    cv.ShowImage("VelY", velY)
    
    cv.WaitKey(0)
