# Workflow:
# import queryContext as context
# [Set whatever variables you want in context]
# context.vars_init()
# context.match(sift, matchdir, [lat], [lon])
# context.characterize()

import os
import os.path
import query
from reader import get_reader

QUERY = None # set before calling characterize()
params = query.PARAMS_DEFAULT.copy()
params.update({
  'checks': 1024,
  'trees': 1,
  'vote_method': 'filter',
  'num_neighbors': 1,
  'confstring': '',
})
count = 0
start = 0
cacheEnable = 0 # instance-local caching of results
ransac_min_filt = 1
print_per = 1
num_images_to_print = 1
put_into_dirs = 0
showHom = 0
locator_function = lambda image: [(image.lat, image.lon)]
cellradius = 236.6
match_callback = None
ambiguity = 75
matchdistance = 25
ncells = 10 # if ambiguity<100, 9 is max possible by geometry
verbosity = 1
resultsdir = os.path.expanduser('~/topmatches')
maindir = os.path.expanduser('/media/DATAPART2')
topnresults = []
initialized = False

# computed based on maindir, QUERY
def vars_init():
    global initialized
    global dbdump
    global dbdir
    global fuzzydir
    global results
    initialized = True
    dbdump = os.path.join(maindir, "Research/collected_images/earthmine-fa10.1,culled/37.871955,-122.270829")
    dbdir = os.path.join(maindir, 'Research/cells/g=100,r=d=236.6/')
    fuzzydir = os.path.join(maindir, 'fuzzylocs/%s' % QUERY)
    results = {}
    for n in topnresults:
        results[n]=0

cache = {}

try:
    if 'NUM_THREADS' in os.environ:
        NUM_THREADS = int(os.environ['NUM_THREADS'])
    else:
        import multiprocessing
        NUM_THREADS = multiprocessing.cpu_count()
    drawtopcorr = 'NO_DRAW' not in os.environ
except:
    import multiprocessing
    NUM_THREADS = multiprocessing.cpu_count()
    drawtopcorr = 1

import shutil
import time

from config import *
from android import AndroidReader
import info
import numpy as np
import corr
import query1GroundTruth
import query2Groundtruth
import groundtruthB
import groundtruthG
import groundtruthO
import groundtruthR
import groundtruthY
import query4GroundTruth
import util

class Img:
    def __init__(self):
        self.lat, self.lon, self.sift = None, None, None

# for test set runs
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

def skew_location(image):
    length = 2*ambiguity
    points = []
    corner = info.moveLocation(image.lat, image.lon, (2**.5)*ambiguity, -45)
    for i in range(length+1):
        row = info.moveLocation(corner[0], corner[1], i, 180)
        for j in range(length+1):
            point = info.moveLocation(row[0],row[1], j, 90)
            if info.distance(image.lat,image.lon, point[0], point[1]) <= ambiguity:
                points.append(point)
    return points

def load_location(image):
    querysift = image.sift
    file = open(os.path.join(fuzzydir, querysift.split("sift.txt")[0])+'.fuz')
    points = []
    for line in file:
        lat, lon = line.strip().split('\t')
        lat = float(lat)
        lon = float(lon)
        points.append((lat,lon))
    return points

def derive_key(closest_cells, querysift):
    return (querysift,) + tuple(sorted(map(lambda (cell, dist): cell, closest_cells)))

# newlat and newlon are skewed locs
def match(siftpath, matchdir, lat, lon, newlat=None, newlon=None):
    assert os.path.basename(siftpath) != siftpath
    assert initialized, "You must call vars_init() first"
    querydir = os.path.dirname(siftpath)
    siftfile = os.path.basename(siftpath)

    # compute closest cells
    closest_cells = util.getclosestcells(newlat or lat, newlon or lon, dbdir)
    cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:ncells] if dist < cellradius + ambiguity + matchdistance]

# Not really needed
#    # query.py filter assumption
#    for cell, dist in cells_in_range:
#        assert cell != '37.8732916946,-122.279128355'

    # cache for fuzz runs
    if cacheEnable:
        key = derive_key(cells_in_range, siftfile)
        if key in cache:
            return cache[key]

    # compute output file paths for the cells
    outputFilePaths = []
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(lat, lon, latcell, loncell)
        outputFilePath = os.path.join(matchdir, siftfile + ',' + cell + ',' + str(actualdist)  + ".res")
        outputFilePaths.append(outputFilePath)

    # start query
    query.run_parallel(dbdir, [c for c,d in cells_in_range], querydir, siftfile, outputFilePaths, params)

    # combine results
    comb_matches = corr.combine_matches(outputFilePaths)
    combined = combine_ransac(comb_matches, ransac_min_filt)

    # top 1
    stats = check_topn_img(siftfile, combined, 1)
    match = any(stats)

    # maybe draw output file
    if drawtopcorr:
        draw_top_corr(querydir, siftfile.split('sift.txt')[0], combined, lat, lon, comb_matches)

    # return statistics and top result
    matchedimg = combined[0][0]
    matches = comb_matches[matchedimg + 'sift.txt']
    if cacheEnable:
        cache[key] = (stats, matchedimg, matches, combined)
    if match_callback:
        match_callback(siftfile, matchdir, stats, matchedimg, combined, cells_in_range)

    # done
    return stats, matchedimg, matches, combined

def draw_top_corr(querydir, query, ranked_matches, qlat, qlon, comb_matches):
    topentry = ranked_matches[0]
    matchedimg = topentry[0]
    score = topentry[1]
    
    dup = "dup" + str(len(ranked_matches) == 1 or score == ranked_matches[1][1])
    
    clat = float(matchedimg.split(",")[0])
    clon = float(matchedimg.split(",")[1][0:-5])
    distance = info.distance(qlat, qlon, clat, clon)

    if put_into_dirs:
        udir = os.path.join(resultsdir, query)
    else:
        udir = resultsdir
    if not os.path.exists(udir):
        os.makedirs(udir)
    queryimgpath = os.path.join(querydir, query + '.jpg')
    if not os.path.exists(queryimgpath):
        queryimgpath = os.path.join(querydir, query + '.JPG')
        assert os.path.exists(queryimgpath)
    i = 0
    for matchedimg, score in ranked_matches[:num_images_to_print]:
        i += 1
        clat = float(matchedimg.split(",")[0])
        clon = float(matchedimg.split(",")[1][0:-5])
        distance = info.distance(qlat, qlon, clat, clon)
        matchimgpath = os.path.join(dbdump, '%s.jpg' % matchedimg)
        match = any(check_img(query + 'sift.txt', ranked_matches[i-1]))
        # rematch for precise fit
        db_matches = comb_matches[matchedimg + 'sift.txt']
        matches = db_matches
        reader = get_reader(params['descriptor'])
        querysiftpath = os.path.join(querydir, query + 'sift.txt')
        matchsiftpath = os.path.join(dbdump, matchedimg + 'sift.txt')
        matches = corr.rematch(reader, querysiftpath, matchsiftpath)
        # concat db matches
        matches.extend(db_matches)

        data = {}
        matchoutpath = os.path.join(udir, query + ';match' + str(i) + ';gt' + str(match)  + ';hom' + str(None) + ';' + matchedimg + '.jpg')
        H, inliers = corr.draw_matches(matches, queryimgpath, matchimgpath, matchoutpath, showHom=showHom, data=data)
        new = os.path.join(udir, query[4:8] + ';match' + str(i) + ';gt' + str(match)  + ';hom' + str(data.get('success')) + ';uniq=' + str(data.get('unique_features')) + ';inliers=' + str(float(sum(inliers))/len(matches)) + ';' + matchedimg + '.jpg')
        os.rename(matchoutpath, new)

        if showHom:
            if put_into_dirs:
                identifier = str(i);
            else:
                identifier = query[4:8] + ':' + str(i)
            H = np.matrix(np.asarray(H))
            with open(os.path.join(udir, 'homography%s.txt' % identifier), 'w') as f:
                print >> f, H
            np.save(os.path.join(udir, 'matches%s.npy' % identifier), matches)

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

def check_img(querysift, entry):
    g,y,r,b,o = 0,0,0,0,0
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
        g += check_truth(querysift.split('sift')[0], entry[0], query4GroundTruth.matches)
    else:
        return [0,0,0,0,0]
    return [g > 0, y > 0, r > 0, b > 0, o > 0]

def check_topn_img(querysift, dupCountLst, topnres=1):
    record = [0]*5
    for entry in dupCountLst[0:topnres]:
        new = check_img(querysift, entry)
        record = map(lambda a,b: a + b, record, new)
    return map(bool, record)

def dump_combined_matches(siftfile, matchdir, stats, matchedimg, matches, cells_in_range):
    # For Aaron's analysis
    table = {}
    for line in open(os.path.join(dbdir, 'cellmap.txt')):
        a, b = line.split()
        table[b] = int(a)
    def cellsetstr(cells):
        cells = sorted(map(lambda (cell, dist): str(table[cell]), cells))
        return '-'.join(cells)
    outputFilePath = os.path.join(matchdir, 'fuzz', siftfile + ',combined,' + cellsetstr(cells_in_range) + ".res")
    d = os.path.dirname(outputFilePath)
    if not os.path.exists(d):
        os.makedirs(d)
    def save(outputFilePath):
        with open(outputFilePath, 'w') as outfile:
            for matchedimg, score in matches:
                outfile.write(str(score))
                outfile.write('\t')
                outfile.write(matchedimg)
                outfile.write('\n')
    save_atomic(save, outputFilePath)
    
def characterize():
    assert QUERY
    assert initialized, "You must call vars_init() first"
    matchdir = os.path.join(maindir, 'Research/results/%s/matchescells(g=100,r=d=236.6),%s,%s' % (QUERY, QUERY, query.searchtype(params)))
    INFO("matchdir=%s" % matchdir)
    querydir = os.path.join(maindir, '%s/' % QUERY)
    global start # XXX
    start = time.time()
    if not os.path.exists(matchdir):
        os.makedirs(matchdir)
    if drawtopcorr:
        if os.path.exists(resultsdir):
            shutil.rmtree(resultsdir)
        os.makedirs(resultsdir)
    reader = make_reader(querydir)
    g_count = 0
    y_count = 0
    r_count = 0
    b_count = 0
    o_count = 0
    global count # XXX
    for image in reader:
        queryfile = image.sift
        for loc in locator_function(image):
            count += 1
            querypath = os.path.join(querydir, queryfile)
            [g, y, r, b, o], matchedimg, matches, combined = match(querypath, matchdir, image.lat, image.lon, loc[0], loc[1])
            # compile statistics
            # top n
            for n in topnresults:
                result = check_topn_img(queryfile, combined, n)
                results[n] += reduce(lambda x,y: x or y, result)
            if count % print_per == 0 and verbosity > 0:
                INFO('speed is %f' % ((time.time()-start)/count))
                for n in topnresults:
                    print "matched {0}\t out of {1}\t in the top {2}\t amb: {3}, ncells:{4}".format(results[n], count, n, ambiguity, ncells)

            if g:
                g_count += 1
                if verbosity > 0 and count % print_per == 0:
                    print "G match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            elif y:
                y_count += 1
                if verbosity > 0 and count % print_per == 0:
                    print "Y match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            elif r:
                r_count += 1
                if verbosity > 0 and count % print_per == 0:
                    print "R match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            elif b:
                b_count += 1
                if verbosity > 0 and count % print_per == 0:
                    print "B match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            elif o:
                o_count += 1
                if verbosity > 0 and count % print_per == 0:
                    print "O match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            else:
                if verbosity > 0 and count % print_per == 0:
                    print "No match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)

    end = time.time()
    elapsed = end - start
    if verbosity > 0:
        print "total time:{0}, avg time:{1}".format(elapsed, elapsed / count)
    total_count = g_count + y_count + r_count + b_count + o_count
    match_rate = float(total_count) / count
    print "g:{0} y:{1} r:{2} b:{3} o:{4} = {5}, out of {6}={7}".format(g_count, y_count, r_count, b_count, o_count, total_count, count, match_rate)

    print "amb: {0}, ncells: {1}".format(ambiguity, ncells)

