import info
import util
import os
import corr
import query
from querySystemCopy import combine_ransac, Img, draw_top_corr

params = query.PARAMS_DEFAULT.copy()
params.update({
  'checks': 1024,
  'trees': 1,
  'vote_method': 'matchonce',
})
cellradius = 236.6
ambiguity = 50
matchdistance = 25
ncells = 7
maindir = os.path.expanduser('~/shiraz')
dbdir = os.path.join(maindir, 'Research/cells/g=100,r=d=236.6/')
matchdir = os.path.expanduser('~/results/%s' % query.searchtype(params))
dbdump = os.path.join(maindir, "Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829")
if not os.path.isdir(matchdir):
    os.makedirs(matchdir)

def match(siftfile, imagefile):
    querydir = os.path.dirname(siftfile)
    siftfile = os.path.basename(siftfile)
    lat, lon = info.getQuerySIFTCoord(siftfile)
    closest_cells = util.getclosestcells(lat, lon, dbdir)
    outputFilePaths = []
    cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:ncells] if dist < cellradius + ambiguity+matchdistance]
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(lat, lon, latcell, loncell)
        outputFilePath = os.path.join(matchdir, siftfile + ',' + cell + ',' + str(actualdist)  + ".res")
        outputFilePaths.append(outputFilePath)
    # start query
    query.run_parallel(dbdir, [c for c,d in cells_in_range], querydir, siftfile, outputFilePaths, params)
    # end query
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(lat, lon, latcell, loncell)
        outputFilePath = os.path.join(matchdir, siftfile + ',' + cell + ',' + str(actualdist)  + ".res")
    comb_matches = corr.combine_matches(outputFilePaths)
    combined = combine_ransac(comb_matches)
    topentry = combined[0]
    matchedimg = topentry[0]
    matches = comb_matches[matchedimg + 'sift.txt']
    return matchedimg, matches

def draw_corr(queryimgpath, matchedimg, matches):
    F, inliers = corr.find_corr(matches)
    matchimgpath = os.path.join(dbdump, '%s.jpg' % matchedimg)
    matchoutpath = os.path.expanduser('~/client-out.png')
    corr.draw_matches(matches, queryimgpath, matchimgpath, matchoutpath, inliers)
    return F, inliers

if __name__ == '__main__':
    sift = os.path.expanduser('~/shiraz/DSC_7595,37.87015,-122.26853sift.txt')
    image = os.path.expanduser('~/shiraz/DSC_7595,37.87015,-122.26853.JPG')
#    sift = os.path.expanduser('~/shiraz/DSC_7638,37.87162,-122.27223sift.txt')
#    image = os.path.expanduser('~/shiraz/DSC_7638,37.87162,-122.27223.JPG')
#    sift = os.path.expanduser('~/shiraz/DSC_7746,37.87203,-122.27022sift.txt')
#    image = os.path.expanduser('~/shiraz/DSC_7746,37.87203,-122.27022.jpg')
#    sift = os.path.expanduser('~/shiraz/DSC_7712,37.87125,-122.26820sift.txt')
#    image = os.path.expanduser('~/shiraz/DSC_7712,37.87125,-122.26820.jpg')
#    sift = os.path.expanduser('~/shiraz/DSC_7765,37.87228,-122.26843sift.txt')
#    image = os.path.expanduser('~/shiraz/DSC_7765,37.87228,-122.26843.JPG')
    matchedimg, matches = match(sift, image)
    draw_corr(image, matchedimg, matches)
