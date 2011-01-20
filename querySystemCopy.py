import os
import os.path
import shutil
import sys
import time

from config import *
from android import AndroidReader
import info
import query
import query1GroundTruth
import groundtruthB
import groundtruthG
import groundtruthO
import groundtruthR
import groundtruthY
import util

QUERY = 'query4'

def parse_result_line(line):
    score = line.split('\t')[0]
    img = line.split('\t')[1].split('sift')[0]
    return score, img

def check_truth(query_str, result_str, groundTruth_dict):
    print query_str
    print result_str
    return result_str in groundTruth_dict[query_str]

def copy_topn_results(imgdir, outdir, filepath, topn=4):
    if  os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    file = open(filepath)
    lines = []
    for line in file:
        lines.append(line)
    for line in lines[0:4]:
        score = line.split('\t')[0]
        img = line.split('\t')[1].split('sift')[0]
        shutil.copyfile(os.path.join(imgdir, ('%s.jpg' % img)), os.path.join(outdir, '%s-%s.jpg' % (score, img)))

def copy_top_match(querydir, query, ranked_matches, match, qlat, qlon):
    topentry = ranked_matches[0]
    matchedimg = topentry[0]
    score = topentry[1]
    
    dup = "dup" + str(len(ranked_matches) == 1 or score == ranked_matches[1][1])
    
    clat = float(matchedimg.split(",")[0])
    clon = float(matchedimg.split(",")[1][0:-5])
    distance = info.distance(qlat, qlon, clat, clon)

    queryimgpath = os.path.join(querydir, query + '.pgm')
    queryoutpath = os.path.join(resultsdir, query + ';query;gt' + str(match)  + ';' + dup + ';' + matchedimg + ';' + str(score) + ';' + str(distance) + '.pgm')
    shutil.copyfile(queryimgpath, queryoutpath)
    for matchedimg, score in ranked_matches:
        if score != topentry[1]:
            break
        clat = float(matchedimg.split(",")[0])
        clon = float(matchedimg.split(",")[1][0:-5])
        distance = info.distance(qlat, qlon, clat, clon)
        matchimgpath = os.path.join(dbdump, '%s.jpg' % matchedimg)
        matchoutpath = os.path.join(resultsdir, query + ';match;gt' + str(match)  + ';' + dup + ';' + matchedimg + ';' + str(score) + ';' + str(distance) + '.jpg')
        shutil.copy(matchimgpath, matchoutpath)

def write_scores(querysift, ranked_matches, outdir):
    if  not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = open(os.path.join(outdir, querysift + ".top"), 'w')
    for matchedimg, score in ranked_matches:
        outfile.write(str(score))
        outfile.write('\t')
        outfile.write(matchedimg)
        outfile.write('\n')

def query2(querydir, querysift, dbdir, mainOutputDir, nClosestCells, copytopmatch, params, lat, lon, copy_top_n_percell=0):
    closest_cells = util.getclosestcells(lat, lon, dbdir)
    outputFilePaths = []
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
        outputFilePath = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(actualdist)  + ".res")
        outputFilePaths.append(outputFilePath)
    # start query
    query.run_parallel(dbdir, [c for c,d in cells_in_range], querydir, querysift, outputFilePaths, params, 8)
    # end query
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(lat, lon, latcell, loncell)
        outputFilePath = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(actualdist)  + ".res")
        if copy_top_n_percell > 0:
            outputDir = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(actualdist))
            copy_topn_results(os.path.join(dbdir, cell), outputDir, outputFilePath, 4)
#    combined = combine_until_dup(outputFilePaths, 1000)
    combined = combine_topn_votes(outputFilePaths, float('inf'))
#    combined = filter_in_range(combined, querysift)
#    write_scores(querysift, combined, "/media/data/combined")
    [g, y, r, b, o] = check_topn_img(querysift, combined, topnresults)
    if copytopmatch:
        match = g or y or r or b or o
        copy_top_match(querydir, querysift.split('sift.txt')[0], combined, match, lat, lon)
    return [g, y, r, b, o]
    
def filter_in_range(ranked_matches, querysift):
    qlat, qlon = info.getQuerySIFTCoord(querysift)
    weighted_matches = []
    for matchedimg, score in ranked_matches:
        clat = float(matchedimg.split(",")[0])
        clon = float(matchedimg.split(",")[1][0:-5])
        distance = info.distance(qlat, qlon, clat, clon)
        if distance < ambiguity+matchdistance:
            weighted_matches.append((matchedimg, 2 * score))
    weighted_matches.sort(key=lambda x: x[1], reverse=True)
    return weighted_matches

def weigh_cells(ranked_matches, dbdir, nClosestCells):
    weighted_matches = []
    for matchedimg, score in ranked_matches:
        lat = float(matchedimg.split(",")[0])
        lon = float(matchedimg.split(",")[1][0:-5])
        closest_cells = util.getclosestcells(lat, lon, dbdir)
        cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:nClosestCells] if dist < cellradius]
        numoverlap = len(cells_in_range)
        if numoverlap == 3:
            weighted_matches.append((matchedimg, score * 4))
        elif numoverlap == 4:
            weighted_matches.append((matchedimg, score * 3))
#        else:
#            print "ERROR! weigh_cells has more overlap than it's supposed to: {0}, {1}".format(numoverlap, matchedimg)
    weighted_matches.sort(key=lambda x: x[1], reverse=True)
    return weighted_matches

def get_top_results(outputFilePath, n):
    file = open(outputFilePath)
    top_results = []
    i = 0
    for line in file:
        if i >= n:
            break
        score, img = parse_result_line(line)
        top_results.append((score, img))
        i += 1
    file.close()
    return top_results

def combine_topn_votes(outputFilePaths, topn):
    #returns true if query if in topn of results
    dupCount = {}
    for outputFilePath in outputFilePaths:
        for score, img in get_top_results(outputFilePath, topn):
            dupCount[img] = dupCount.get(img, 0) + float(score)
    dupCountLst = dupCount.items()
#    dupCountLst.sort(key=lambda x: x[1])
#    dupCountLst.reverse()
    dupCountLst.sort(key=lambda x: x[1], reverse=True)
    return dupCountLst

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
        elif QUERY == 'query4':
            pass
        else:
            assert False

    return [g > 0, y > 0, r > 0, b > 0, o > 0]
    
def skew_location(querysift, radius):
    center = info.getQuerySIFTCoord(querysift)
    length = 2*radius
    points = []
    corner = info.moveLocation(center[0], center[1], (2**.5)*radius, -45)
    for i in range(length+1):
        row = info.moveLocation(corner[0], corner[1], i, 180)
        for j in range(length+1):
            point = info.moveLocation(row[0],row[1], j, 90)
            if info.distance(center[0],center[1], point[0], point[1]) <= radius:
                points.append(point)
    return points
    
def characterize(querydir, dbdir, mainOutputDir, n, copytopmatch, params):
    start = time.time()
    if not os.path.exists(mainOutputDir):
        os.makedirs(mainOutputDir)
    if copytopmatch:
        if os.path.exists(resultsdir):
            shutil.rmtree(resultsdir)
        os.makedirs(resultsdir)
    reader = AndroidReader(querydir)
    g_count = 0
    y_count = 0
    r_count = 0
    b_count = 0
    o_count = 0
    count = 0
    for image in reader:
        queryfile = image.sift
        count += 1
        [g, y, r, b, o] = query2(querydir, queryfile, dbdir, mainOutputDir, n, copytopmatch, params, image.lat, image.lon)
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

def get_num_imgs_in_range(range, celldir, querydir='/home/zhangz/.gvfs/data on 128.32.43.40/query3/'):
    queries = util.getSiftFileNames(querydir)
    total_in_range=0
    for querysift in queries:
        lat, lon = info.getQuerySIFTCoord(querysift)
        cell, dist = util.getclosestcell(lat, lon, celldir)
        numinrange = util.getNumJPGInRange(lat, lon, os.path.join(celldir,cell), range)
        total_in_range+=numinrange
    print "average # of imgs within {0} meters of query: {1}".format(range, float(total_in_range)/len(queries))

cellradius = 236.6
ambiguity = 50
matchdistance = 25
ncells = 7   #if ambiguity<100, 7 is max possible by geometry
topnresults = 1
verbosity = 1
copytopmatch = True
resultsdir = os.path.expanduser('~/topmatches')
maindir = os.path.expanduser('~/.gvfs/data on 128.32.43.40')
params = query.PARAMS_DEFAULT.copy()
params.update({
  'checks': 512,
  'trees': 1,
  'distance_type': 'euclidean',
  'vote_method': 'highest',
  'confstring': '',
})
dbdump = os.path.join(maindir, "Research/collected_images/earthmine-new,culled/37.871955,-122.270829")
if __name__ == "__main__":
#    querydir = os.path.join(maindir, 'Research/collected_images/query/%s/' % QUERY)
#    querydir = os.path.join(maindir, '%s/' % QUERY)
    querydir = os.path.join(maindir, 'Research/collected_images/droid/%s/' % QUERY)
    dbdir = os.path.join(maindir, 'Research/cellsg=100,r=d=236.6/')
    matchdir = os.path.join(maindir, 'Research/results(%s)/matchescells(g=100,r=d=236.6),%s,%s' % (QUERY, QUERY, query.searchtype(params)))
    if len(sys.argv) > 4:
        print "USAGE: {0} QUERYDIR DBDIR OUTPUTDIR".format(sys.argv[0])
        sys.exit()
    elif len(sys.argv) == 4:
        querydir = sys.argv[1]
        dbdir = sys.argv[2]
        matchdir = sys.argv[3]
    topnresults = 1
    INFO("matchdir=%s" % matchdir)
    characterize(querydir, dbdir, matchdir, ncells, copytopmatch, params)
