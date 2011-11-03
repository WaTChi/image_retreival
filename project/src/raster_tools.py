from PIL import Image
import colorsys

class RasterImage(object):
    def __init__(self, img):
        self.im = Image.open(img).load()

    def sampleHSV(self, x, y, size):
        """returns (h,s,v) each of type uint8"""
        accum = []
        # swap
        y, x = int(x), int(y)
        radius = (size - 1)/2

        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                x2, y2 = x + dx, y + dy
                try:
                    accum.append(self.im[x2, y2])
                except IndexError:
                    pass
        mean = reduce(
            lambda a,b: (a[0] + b[0], a[1] + b[1], a[2] + b[2]),
            accum)
        a = len(accum)*255.0
        mean = (mean[0]/a, mean[1]/a, mean[2]/a)
        hsv = colorsys.rgb_to_hsv(*mean)
        p = (int(hsv[0]*255),
             int(hsv[1]*255),
             int(hsv[2]*255))
        return p

    def getDepth(self, x, y):
        """returns float32 (meters)"""
        raise NotImplementedError

    def getPlane(self, x, y):
        """returns uint32"""
        raise NotImplementedError
