from context import _Context

class _SingleImageContext(_Context):
    
    def __init__(self, context=None):
        super(_SingleImageContext, self).__init__(context)
        self.trackingdir = None;

    @property
    def querydir(self):
        raise Exception("No querydir for SingleImageContext")

    def iter_queries_unfiltered(self):
        def iter():
            firstfile = sorted([os.path.join(trackingdir, file) for file in os.listdir(trackingdir)
                                if os.path.splitext(file)[1] in (".jpg", ".JPG")])[0]
            image = _Query()
            image.jpgpath = os.path.join(trackingdir, firstfile)
            image.siftpath = os.path.join(trackingdir, os.path.splitext(file)[0] + "sift.txt")
            lat, lon = self.read_coords()
            image.setSensorCoord(lat, lon)
            image.check()
            yield image
        return iter()

    def read_coords(self):
        coords_file = os.path.join(trackingdir, "coords.txt")
        if os.path.exists(coords_file):
            lines = open(coords_file).readlines()
            return float(lines[1]), float(lines[2])


DEFAULT_CONTEXT = _SingleImageContext().frozen_copy()
