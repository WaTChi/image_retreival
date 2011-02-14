import cv
import time

MAX_CORNERS = 500;

if __name__ == '__main__':
    start = time.time()
    imgA = cv.LoadImage("/home/etzeng/tracking/imga.jpg",cv.CV_LOAD_IMAGE_GRAYSCALE)
    imgB = cv.LoadImage("/home/etzeng/tracking/imgb.jpg",cv.CV_LOAD_IMAGE_GRAYSCALE)
    
    img_sz = cv.GetSize(imgA)
    win_size = 10
    
    imgC = cv.LoadImage("/home/etzeng/tracking/test/OpticalFlowB.jpg",cv.CV_LOAD_IMAGE_UNCHANGED)
    
    eig_image = cv.CreateImage(img_sz, cv.IPL_DEPTH_32F, 1)
    tmp_image = cv.CreateImage(img_sz, cv.IPL_DEPTH_32F, 1)
    
    corner_count = MAX_CORNERS

#  CvPoint2D32f* cornersA     = new CvPoint2D32f[ MAX_CORNERS ];
  
    cornersA = cv.GoodFeaturesToTrack(imgA, eig_image, tmp_image, corner_count, 0.01, 5.0, None, 3, 0, 0.04)
    cornersA = cv.FindCornerSubPix(imgA, cornersA, (win_size, win_size), (-1, -1), (cv.CV_TERMCRIT_ITER | cv.CV_TERMCRIT_EPS, 20, 0.03))
  
    # Call the Lucas Kanade algorithm
    features_found = []
    feature_errors = []

    pyr_sz = (imgA.width + 8, imgB.height / 3)
    
    pyrA = cv.CreateImage(pyr_sz, cv.IPL_DEPTH_32F, 1)
    pyrB = cv.CreateImage(pyr_sz, cv.IPL_DEPTH_32F, 1)
    
    cornersB = []
#  CvPoint2D32f* cornersB = new CvPoint2D32f[ MAX_CORNERS ];
    (cornersB, features_found, feature_errors) = cv.CalcOpticalFlowPyrLK(imgA, imgB, pyrA, pyrB, cornersA, (win_size, win_size), 5, (cv.CV_TERMCRIT_ITER | cv.CV_TERMCRIT_EPS, 20, .3), 0)
  
    print "TIME ELAPSED: ", str(time.time() - start)

    for i in range(corner_count):
        if features_found[i] == 0 or feature_errors[i] > 550:
            print "Error is {0}".format(feature_errors[i])
            continue
        print "Got it"
        p0 = (cv.Round(cornersA[i][0]), cv.Round(cornersA[i][1]))
        p1 = (cv.Round(cornersB[i][0]), cv.Round(cornersB[i][1]))
        cv.Line(imgC, p0, p1, cv.CV_RGB(255,0,0), 2)

    cv.NamedWindow("ImageA", 0)
    cv.NamedWindow("ImageB", 0)
    cv.NamedWindow("LKpyr_OpticalFlow",0)
    
    cv.ShowImage("ImageA", imgA)
    cv.ShowImage("ImageB", imgB)
    cv.ShowImage("LKpyr_OpticalFlow", imgC)
    
    cv.WaitKey(0)
