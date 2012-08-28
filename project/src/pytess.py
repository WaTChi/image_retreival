import cv2.cv as cv
import tesseract
import sys
import os


def separatePath(path):
    splitat='/'
    myList = ''.join([ s if s not in splitat else ' ' for s in path]).split()
    return myList

def recognize(image):
    api = tesseract.TessBaseAPI()
    api.Init(".","eng",tesseract.OEM_DEFAULT)
    api.SetPageSegMode(tesseract.PSM_SINGLE_BLOCK)
    tesseract.SetCvImage(image,api)
    text=api.GetUTF8Text()
    conf=api.MeanTextConf()
    image=None
    print text
    print "Confidence: " + str(conf)
    return text, conf

if __name__ == '__main__':
    mode = int(sys.argv[2])
    
    if mode == 1:
#        myList = separatePath(sys.argv[1])
        os.system("matlab -nodisplay -nodesktop -nojvm -nosplash -r \"correct_image('', '" + sys.argv[1] + "'); exit\"")
#        image = myList[:-1] + "corrected" + myList[-1]
        image = "corrected"+sys.argv[1] 
        image = cv.LoadImage(image, cv.CV_LOAD_IMAGE_GRAYSCALE)
        t2, c2 = recognize(image)
    
    image = cv.LoadImage(sys.argv[1], cv.CV_LOAD_IMAGE_GRAYSCALE)
    t, c = recognize(image)
    
    if mode == 1 and c2 > c:
        t = t2
        c = c2
    
    f = open(sys.argv[1] + '.txt', 'wb')
    f.write(t + '##' + str(c) + '\n')
    f.close()

#f = open('answers.txt', 'wb')
#
#for x in xrange(250):
#    f.write(str(x+1) + ',\n')
#f.close()