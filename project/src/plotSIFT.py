import sys
from PIL import Image, ImageDraw
import math
import os
import numpy as np

def plotvisualwords(siftfile, imagefile, outfile):
    
    feats, coords, size = readfile(siftfile)


    #Draw on old image, display, and save:
    im = Image.open(imagefile)
    draw = ImageDraw.Draw(im)
    for i in range(len(coords)):
        draw.ellipse(((coords[i][0]-size[i],coords[i][1]-size[i]),
                    (coords[i][0]+size[i],coords[i][1]+size[i])),
                    outline="hsl(20,100%,50%)")

    im.save(outfile)
    
def readfile(path):
    
    if os.path.exists(path):
        f = open(path, 'r')
        (numfeats, dim) = f.readline().split()
        numfeats = int(numfeats)
        dim = int(dim)
        feats = []
        coords = []
        sizes = []
        for i in range(numfeats):
            #Read first line...scale, orientation, etc
            l = f.readline().strip().split()
            y = l[0]
            x = l[1]
            s = float(l[2])
            coords.append((int(round(float(x))),int(round(float(y)))))
            sizes.append(s)
            count = 0
            desc = []
            while count < dim:
                l = f.readline().strip().split()
                count += len(l)
                desc.extend(l)
            try:
                feats.append(np.array(desc, dtype=np.uint8))
            except NameError:
                pass
        feats = np.array(feats, dtype=np.uint8)
        f.close()
        return feats, coords, sizes
    else:
        raise OSError("{p} does not exist.".format(p=path))


    return feats, coords, sizes

if __name__=="__main__":
    plotvisualwords(sys.argv[1], sys.argv[2], sys.argv[3])
