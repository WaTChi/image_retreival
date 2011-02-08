import os
import os.path
import shutil
import sys
import time

from config import *
from android import AndroidReader
import info
import numpy as np
import query
import corr
import query1GroundTruth
import query2Groundtruth
import groundtruthB
import groundtruthG
import groundtruthO
import groundtruthR
import groundtruthY
import util

QUERY = 'query3'
try:
    if 'NUM_THREADS' in os.environ:
        NUM_THREADS = int(os.environ['NUM_THREADS'])
    else:
        import multiprocessing
        NUM_THREADS = multiprocessing.cpu_count()
    drawtopcorr = 'NO_DRAW' not in os.environ
    drawfailures = 'DRAW_FAIL' in os.environ
except:
    import multiprocessing
    NUM_THREADS = multiprocessing.cpu_count()
    drawtopcorr = 1
    drawfailures = 0

class Img:
    def __init__(self):
        self.lat, self.lon, self.sift = None, None, None

def make_reader(querydir):
    if QUERY == 'query4':
        return AndroidReader(querydir)
    INFO(querydir)
    def iter():
        for file in util.getSiftFileNames(querydir):
            image = Img()
            image.sift = file
            image.lat, image.lon = info.getQuerySIFTCoord(file)
            yield image
    return iter()

def check_truth(query_str, result_str, groundTruth_dict):
    return result_str in groundTruth_dict[query_str]

def draw_top_corr(querydir, query, ranked_matches, match, qlat, qlon, comb_matches):
    topentry = ranked_matches[0]
    matchedimg = topentry[0]
    score = topentry[1]
    
    dup = "dup" + str(len(ranked_matches) == 1 or score == ranked_matches[1][1])
    
    clat = float(matchedimg.split(",")[0])
    clon = float(matchedimg.split(",")[1][0:-5])
    distance = info.distance(qlat, qlon, clat, clon)

    udir = os.path.join(resultsdir, str(match))
#    udir = os.path.join(resultsdir, str(match), query)
    if not os.path.exists(udir):
        os.makedirs(udir)
    queryimgpath = os.path.join(querydir, query + '.jpg')
    i = 0
    for matchedimg, score in ranked_matches[:1]:
#        if score != topentry[1]:
#            break
#        i += 1 # XXX tmp image analysis code to spit out top 10
#        if i > 9:
#            break;
        clat = float(matchedimg.split(",")[0])
        clon = float(matchedimg.split(",")[1][0:-5])
        distance = info.distance(qlat, qlon, clat, clon)
        matchimgpath = os.path.join(dbdump, '%s.jpg' % matchedimg)
        matches = comb_matches[matchedimg + 'sift.txt']
        F, inliers = corr.find_corr(matches)
        matchoutpath = os.path.join(udir, query + ';match' + str(i) + '(' + str(int(score)) + ');gt' + str(match)  + ';' + dup + ';' + matchedimg + ';' + str(score) + ';' + str(distance) + '.jpg')
        corr.draw_matches(matches, queryimgpath, matchimgpath, matchoutpath, inliers)

def query2(querydir, querysift, dbdir, mainOutputDir, nClosestCells, drawtopcorr, params, lat, lon):
    closest_cells = util.getclosestcells(lat, lon, dbdir)
    outputFilePaths = []
    print len(closest_cells)
    cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:nClosestCells] if dist < cellradius + ambiguity+matchdistance]
    if verbosity > 0:
        print "checking query: {0} \t against {1} \t cells".format(querysift, len(cells_in_range))
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(lat, lon, latcell, loncell)
        if verbosity > 1:
            print "querying cell: {0}, distance: {1} with:{2}".format(cell, actualdist, querysift)
        print querysift + ',' + cell + ',' + str(actualdist)  + ".res"
        outputFilePath = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(actualdist)  + ".res")
        outputFilePaths.append(outputFilePath)
    # start query
    query.run_parallel(dbdir, [c for c,d in cells_in_range], querydir, querysift, outputFilePaths, params, NUM_THREADS)
    # end query
    comb_matches = corr.combine_matches(outputFilePaths)
    combined = combine_ransac(comb_matches)
    [g, y, r, b, o] = check_topn_img(querysift, combined, topnresults)
    match = g or y or r or b or o
    if drawtopcorr or (not match and drawfailures):
        draw_top_corr(querydir, querysift.split('sift.txt')[0], combined, match, lat, lon, comb_matches)
    return [g, y, r, b, o]

def combine_ransac(counts):
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    filtered = {}
    bound = -1
    num_filt = 0
    for siftfile, matches in sorted_counts:
      siftfile = siftfile[:-8]
      if len(matches) < bound or num_filt > 20:
        INFO('stopped after filtering %d' % num_filt)
        break
      num_filt += 1
      F, inliers = corr.find_corr(matches)
      bound = max(sum(inliers), bound)
      pts = np.ndarray(len(matches), np.object)
      pts[0:len(matches)] = matches
      if any(inliers):
        filtered[siftfile] = list(np.compress(inliers, pts))
    rsorted_counts = sorted(filtered.iteritems(), key=lambda x: len(x[1]), reverse=True)
    def condense(list):
        return map(lambda x: (x[0], len(x[1])), list)
    def condense2(list):
        return map(lambda x: (x[0][:-8], len(x[1])), list)
    if not rsorted_counts:
      INFO('W: postcomb ransac rejected everything, not filtering')
      return condense2(sorted_counts)
    return condense(rsorted_counts)

def check_topn_img(querysift, dupCountLst, topnres=1):
    g = 0
    y = 0
    r = 0
    b = 0
    o = 0
    for entry in dupCountLst[0:topnres]:
        if QUERY == 'query1':
            g += check_truth(querysift.split('sift')[0], entry[0], query1GroundTruth.matches)
        elif QUERY == 'query3' or QUERY == 'queryeric':
            g += check_truth(querysift.split('sift')[0], entry[0], groundtruthG.matches)
            y += check_truth(querysift.split('sift')[0], entry[0], groundtruthY.matches)
            r += check_truth(querysift.split('sift')[0], entry[0], groundtruthR.matches)
            b += check_truth(querysift.split('sift')[0], entry[0], groundtruthB.matches)
            o += check_truth(querysift.split('sift')[0], entry[0], groundtruthO.matches)
        elif QUERY == 'query2':
            g += check_truth(querysift.split('sift')[0], entry[0], query2Groundtruth.matches)
        elif QUERY == 'query4':
            pass
        else:
            assert False
    return [g > 0, y > 0, r > 0, b > 0, o > 0]
    
def characterize(querydir, dbdir, mainOutputDir, n, drawtopcorr, params):
    start = time.time()
    if not os.path.exists(mainOutputDir):
        os.makedirs(mainOutputDir)
    if drawtopcorr or drawfailures:
        if os.path.exists(resultsdir):
            shutil.rmtree(resultsdir)
        os.makedirs(resultsdir)
    reader = make_reader(querydir)
    g_count = 0
    y_count = 0
    r_count = 0
    b_count = 0
    o_count = 0
    count = 0
    for image in reader:
        queryfile = image.sift
        count += 1
        [g, y, r, b, o] = query2(querydir, queryfile, dbdir, mainOutputDir, n, drawtopcorr, params, image.lat, image.lon)
        if g:
            g_count += 1
            if verbosity > 0:
                print "G match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        elif y:
            y_count += 1
            if verbosity > 0:
                print "Y match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        elif r:
            r_count += 1
            if verbosity > 0:
                print "R match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        elif b:
            b_count += 1
            if verbosity > 0:
                print "B match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        elif o:
            o_count += 1
            if verbosity > 0:
                print "O match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        else:
            if verbosity > 0:
                print "No match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
    end = time.time()
    elapsed = end - start
    if verbosity > 0:
        print "total time:{0}, avg time:{1}".format(elapsed, elapsed / count)
    total_count = g_count + y_count + r_count + b_count + o_count
    match_rate = float(total_count) / count
    print "g:{0} y:{1} r:{2} b:{3} o:{4} = {5}, out of {6}={7}".format(g_count, y_count, r_count, b_count, o_count, total_count, count, match_rate)

cellradius = 236.6
ambiguity = 75
matchdistance = 25
ncells = 8   #if ambiguity<100, 8 is max possible by geometry
topnresults = 1
verbosity = 1
resultsdir = os.path.expanduser('~/shiraz/topmatches')
maindir = os.path.expanduser('~/shiraz')
params = query.PARAMS_DEFAULT.copy()
params.update({
  'checks': 1024,
  'trees': 1,
  'distance_type': 'euclidean',
  'vote_method': 'filter',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})
dbdump = os.path.join(maindir, "Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829")
if __name__ == "__main__":
    querydir = os.path.join(maindir, '%s/' % QUERY)
    dbdir = os.path.join(maindir, 'Research/cells/g=100,r=d=236.6/')
    matchdir = os.path.join(maindir, 'Research/results/%s/matchescells(g=100,r=d=236.6),%s,%s' % (QUERY, QUERY, query.searchtype(params)))
    if len(sys.argv) > 4:
        print "USAGE: {0} QUERYDIR DBDIR OUTPUTDIR".format(sys.argv[0])
        sys.exit()
    elif len(sys.argv) == 4:
        querydir = sys.argv[1]
        dbdir = sys.argv[2]
        matchdir = sys.argv[3]
    topnresults = 1
    INFO("matchdir=%s" % matchdir)
    characterize(querydir, dbdir, matchdir, ncells, False, params)
