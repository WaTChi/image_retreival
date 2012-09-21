import util

#inputDir = "/media/DATAPART2/Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829"
#outputDir = "/media/DATAPART1/earthmine-fa10.1-culled,r=%s,d=%s/" % (str(r), str(d))

for x in [150,350]:
    d = r = x
    inputDir = "/media/DATAPART1/oakland/earthmine/rect"
    outputDir = "/media/DATAPART1/oak,r=%s,d=%s/" % (str(r), str(d))

    util.makecells(lat=37.813, lon=-122.28, length=2500, inputDir=inputDir, distance=d, radius=r, outputDir=outputDir)
