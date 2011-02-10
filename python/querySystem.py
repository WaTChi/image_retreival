import os
import os.path
import shutil
import sys
import time

HOME = os.path.expanduser('~')
QUERY = 'query2'
ncells = 8

from config import *
import groundtruthB
import groundtruthG
import groundtruthO
import groundtruthR
import groundtruthY
import info
import query1GroundTruth
import query2Groundtruth
import util
import query
import datetime
import corr
import numpy as np

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

def parse_result_line(line):
    score = line.split('\t')[0]
    img = line.split('\t')[1].split('sift')[0]
    return score, img

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

def parse_result_line(line):
    score = line.split('\t')[0]
    img = line.split('\t')[1].split('sift')[0]
    return score, img

def check_truth(query_str, result_str, groundTruth_dict):
    return result_str in groundTruth_dict[query_str]

def derive_key(closest_cells, querysift):
    return (querysift,) + tuple(sorted(map(lambda (cell, dist): cell, closest_cells)))

cache = {}
def run_query(newlat, newlon, querydir, querysift, dbdir, mainOutputDir, nClosestCells, copytopmatch, closest_cells, params, copy_top_n_percell=0):
    cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:nClosestCells] if dist < cellradius + ambiguity+matchdistance]
    key = derive_key(cells_in_range, querysift)
    if key in cache:
        return cache[key]
    outputFilePaths = []
#    # query.py filter assumption
#    # I think it doesn't matter too much since ties are
#    # resolved in favor of the query feature.
#    for cell, dist in cells_in_range:
#        assert cell != '37.8732916946,-122.279128355'
    latquery, lonquery = info.getQuerySIFTCoord(querysift)
    if verbosity > 0:
        print "checking query: {0} \t against {1} \t cells".format(querysift, len(cells_in_range))
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(latquery, lonquery, latcell, loncell)
        if verbosity > 1:
            print "querying cell: {0}, distance: {1} with:{2}".format(cell, actualdist, querysift)
        outputFilePath = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(actualdist)  + ".res")
        outputFilePaths.append(outputFilePath)
#     start query
    query.run_parallel(dbdir, [c for c,d in cells_in_range], querydir, querysift, outputFilePaths, params, NUM_THREADS)
#     end query
    comb_matches = corr.combine_matches(outputFilePaths)
    combined = combine_ransac(comb_matches, 9999999)

    # For Aaron's analysis
    table = {}
    for line in open(os.path.join(dbdir, 'cellmap.txt')):
        a, b = line.split()
        table[b] = int(a)
    def cellsetstr(cells):
        cells = sorted(map(lambda (cell, dist): str(table[cell]), cells))
        return '-'.join(cells)
    outputFilePath = os.path.join(mainOutputDir, 'fuzz', querysift + ',combined,' + cellsetstr(cells_in_range) + ".res")
    d = os.path.dirname(outputFilePath)
    if not os.path.exists(d):
        os.makedirs(d)
    def save(outputFilePath):
        with open(outputFilePath, 'w') as outfile:
            for matchedimg, score in combined:
                outfile.write(str(score))
                outfile.write('\t')
                outfile.write(matchedimg)
                outfile.write('\n')
    save_atomic(save, outputFilePath)

    results = {}
    for n in topnresults:
        result = check_topn_img(querysift, combined, n)
        results[n] = reduce(lambda x,y: x or y, result)
    cache[key] = results
    return results

def combine_ransac(counts, min_filt=0):
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    filtered = {}
    bound = -1
    num_filt = 0
    for siftfile, matches in sorted_counts:
      siftfile = siftfile[:-8]
      if num_filt > min_filt and (len(matches) < bound or num_filt > 20):
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

def skew_location(querysift, radius):
    center = info.getQuerySIFTCoord(querysift)
    return _skew_location(center, radius)

def _skew_location(center, radius):
    length = 2*radius
    points = []
    corner = info.moveLocation(center[0], center[1], (2**.5)*radius, -45)
    for i in range(length+1):
        row = info.moveLocation(corner[0], corner[1], i, 180)
        for j in range(length+1):
            point = info.moveLocation(row[0],row[1], j, 90)
            if info.distance(center[0],center[1], point[0], point[1]) <= radius:
                #newquerysift = querysift.split(',')[0]+','+str(point[0])+','+str(point[1])+'sift.txt'
                points.append(point)
    return points

def load_locations(querysift):
    file = open(os.path.join(fuzzydir, querysift.split("sift.txt")[0])+'.fuz')
    points = []
    for line in file:
        lat, lon = line.strip().split('\t')
        lat = float(lat)
        lon = float(lon)
        points.append((lat,lon))
    return points

def characterize_fuzzy(querydir, dbdir, mainOutputDir, ncells, params):
    start = time.time()
    if not os.path.exists(mainOutputDir):
        os.makedirs(mainOutputDir)
    files = util.getSiftFileNames(querydir)
    results={}
    for n in topnresults:
        results[n]=0
    count = 0
    start = time.time()
    for queryfile in files:
        print "checking: {0}".format(queryfile)
#        for newlat, newlon in [info.getQuerySIFTCoord(queryfile)]:
#        for newlat, newlon in skew_location(queryfile, ambiguity):
        for newlat, newlon in load_locations(queryfile):
            closest_cells = util.getclosestcells(newlat, newlon, dbdir)
            count += 1
            result = run_query(newlat, newlon, querydir, queryfile, dbdir, mainOutputDir, ncells, copytopmatch, closest_cells, params)
            for n in topnresults:
                results[n]+=result[n]
            if verbosity > 0:
                INFO('speed is %f' % ((time.time()-start)/count))
                for n in topnresults:
                    print "matched {0}\t out of {1}\t in the top {2}\t amb: {3}, ncells:{4}".format(results[n], count, n, ambiguity, ncells)
        print datetime.datetime.today().strftime("time: %l:%M:%S")
        for n in topnresults:
            print "matched {0}\t out of {1}\t in the top {2}\t amb: {3}, ncells:{4}".format(results[n], count, n, ambiguity, ncells)

    end = time.time()
    elapsed = end - start
    if verbosity > 0:
        print "total time:{0}, avg time:{1}".format(elapsed, elapsed / count)

    print "amb: {0}, ncells: {1}".format(ambiguity, ncells)

    for n in topnresults:
        total_count = results[n]
        match_rate = float(total_count) / count
        print "matched {0}\t out of {1}\t = {2}\t in the top {3}".format(total_count, count, match_rate, n)

cellradius = 236.6
ambiguity = 75
matchdistance = 25
topnresults = [1,2,5,10]
verbosity = 0
copytopmatch = False
resultsdir = '/media/data/topmatches'
maindir = HOME + "/shiraz"
fuzzydir = os.path.join(maindir, 'fuzzylocs/%s' % QUERY)
dbdump = os.path.join(maindir, "Research/collected_images/earthmine-new,culled/37.871955,-122.270829")
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
if __name__ == "__main__":
    querydir = os.path.join(maindir, '%s/' % QUERY)
    dbdir = os.path.join(maindir, 'Research/cells/g=100,r=d=236.6/')
   #for gorgan
    matchdir = os.path.join(maindir, 'Research/results/%s/matchescells(g=100,r=d=236.6),%s,%s' % (QUERY, QUERY, query.searchtype(params)))
    if len(sys.argv) > 4:
        print "USAGE: {0} QUERYDIR DBDIR OUTPUTDIR".format(sys.argv[0])
        sys.exit()
    elif len(sys.argv) == 4:
        querydir = sys.argv[1]
        dbdir = sys.argv[2]
        matchdir = sys.argv[3]
    for n in [75]:
        ambiguity = n
        print "ambiguity:{0}".format(n)
        print querydir
        print dbdir
        print matchdir
        characterize_fuzzy(querydir, dbdir, matchdir, ncells, params)
        print "ambiguity:{0}".format(n)
        print querydir
        print dbdir
        print matchdir
