import android

BASEDIR = "/media/DATAPART2/query5horizontal"
THRESHOLD = 5

if __name__ == '__main__':
    ar = android.AndroidReader(BASEDIR)
    print [img.jpg[:-4] for img in ar if img.gps_acc <= 5]
    #for i in sorted([img.gps_acc for img in ar if img.gps_acc <= 5]):
    #    print i
