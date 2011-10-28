import os

class AbstractColorMap(object):
    def __init__(self, image_path):
        self.image_path = image_path

    def getColor(self, x, y):
        raise NotImplementedError

class ColorMap(AbstractColorMap):
    def getRGB(self, x, y):
        """returns (r, g, b) at point"""
        pass # TODO

    def getHSV(self, x, y):
        """returns (h, s, v) at point"""
        pass # TODO

class PlaneMap(AbstractColorMap):
    def getPlane(self, x, y):
        """returns unique identifier for a plane"""
        pass # TODO

class DepthMap(AbstractColorMap):
    def __init__(self, image_path, info_info):
        AbstractColorMap.__init__(self, image_path)
        self.image_info = image_info

    def getCoord(self, x, y):
        """returns (lat, lon, alt)"""
        pass # TODO

    def getDepth(self, x, y):
        """returns distance from camera in meters"""
        pass # TODO

