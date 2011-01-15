import os
import os.path
import shutil
import sys
import time

from config import *
import info
import query
import query1GroundTruth
import query2Groundtruth
import groundtruthB
import groundtruthG
import groundtruthO
import groundtruthR
import groundtruthY
import util

QUERY = 'query3'

def parse_result_line(line):
    score = line.split('\t')[0]
    img = line.split('\t')[1].split('chog')[0]
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
        img = line.split('\t')[1].split('chog')[0]
        shutil.copyfile(os.path.join(imgdir, ('%s.jpg' % img)), os.path.join(outdir, '%s-%s.jpg' % (score, img)))

def copy_top_match(querydir, query, ranked_matches, match):
    topentry = ranked_matches[0]
    matchedimg = topentry[0]
    score = topentry[1]
    
    dup = "dup" + str(len(ranked_matches) == 1 or score == ranked_matches[1][1])
    
    qlat = float(query.split(',')[1])
    qlon = float(query.split(',')[2])
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

def write_scores(querychog, ranked_matches, outdir):
    if  not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = open(os.path.join(outdir, querychog + ".top"), 'w')
    for matchedimg, score in ranked_matches:
        outfile.write(str(score))
        outfile.write('\t')
        outfile.write(matchedimg)
        outfile.write('\n')

def query2(querydir, querychog, dbdir, mainOutputDir, nClosestCells, copytopmatch, params, copy_top_n_percell=0):
    lat, lon = info.getQueryCHOGCoord(querychog)
    closest_cells = util.getclosestcells(lat, lon, dbdir)
    outputFilePaths = []
    cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:nClosestCells] if dist < cellradius + ambiguity+matchdistance]
    latquery, lonquery = info.getQueryCHOGCoord(querychog)
    if verbosity > 0:
        print "checking query: {0} \t against {1} \t cells".format(querychog, len(cells_in_range))
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(latquery, lonquery, latcell, loncell)
        if verbosity > 1:
            print "querying cell: {0}, distance: {1} with:{2}".format(cell, actualdist, querychog)
        outputFilePath = os.path.join(mainOutputDir, querychog + ',' + cell + ',' + str(actualdist)  + ".res")
        outputFilePaths.append(outputFilePath)
    # start query
    query.run_parallel(dbdir, [c for c,d in cells_in_range], querydir, querychog, outputFilePaths, params)
    # end query
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(latquery, lonquery, latcell, loncell)
        outputFilePath = os.path.join(mainOutputDir, querychog + ',' + cell + ',' + str(actualdist)  + ".res")
        if copy_top_n_percell > 0:
            outputDir = os.path.join(mainOutputDir, querychog + ',' + cell + ',' + str(actualdist))
            copy_topn_results(os.path.join(dbdir, cell), outputDir, outputFilePath, 4)
#    combined = combine_until_dup(outputFilePaths, 1000)
    combined = combine_topn_votes(outputFilePaths, float('inf'))
#    combined = filter_in_range(combined, querychog)
#    write_scores(querychog, combined, "/media/data/combined")
    [g, y, r, b, o] = check_topn_img(querychog, combined, topnresults)
    if copytopmatch:
        match = g or y or r or b or o
        copy_top_match(querydir, querychog.split('chog.txt')[0], combined, match)
    return [g, y, r, b, o]
    
def filter_in_range(ranked_matches, querychog):
    qlat, qlon = info.getQueryCHOGCoord(querychog)
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

def check_topn_img(querychog, dupCountLst, topnres=1):
    g = 0
    y = 0
    r = 0
    b = 0
    o = 0
    for entry in dupCountLst[0:topnres]:
        if QUERY == 'query1':
            g += check_truth(querychog.split('chog')[0], entry[0], query1GroundTruth.matches)
        elif QUERY == 'query2':
            g += check_truth(querychog.split('chog')[0], entry[0], query2Groundtruth.matches)
        elif QUERY == 'query3':
            g += check_truth(querychog.split('chog')[0], entry[0], groundtruthG.matches)
            y += check_truth(querychog.split('chog')[0], entry[0], groundtruthY.matches)
            r += check_truth(querychog.split('chog')[0], entry[0], groundtruthR.matches)
            b += check_truth(querychog.split('chog')[0], entry[0], groundtruthB.matches)
            o += check_truth(querychog.split('chog')[0], entry[0], groundtruthO.matches)
        else:
            assert False

    return [g > 0, y > 0, r > 0, b > 0, o > 0]
    
def characterize(querydir, dbdir, mainOutputDir, n, copytopmatch, params):
    start = time.time()
    if not os.path.exists(mainOutputDir):
        os.makedirs(mainOutputDir)
    if copytopmatch:
        if os.path.exists(resultsdir):
            shutil.rmtree(resultsdir)
        os.makedirs(resultsdir)
    files = util.getCHOGFileNames(querydir)
    g_count = 0
    y_count = 0
    r_count = 0
    b_count = 0
    o_count = 0
    count = 0
    for queryfile in files:
            count += 1
            [g, y, r, b, o] = query2(querydir, queryfile, dbdir, mainOutputDir, n, copytopmatch, params)
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
    queries = util.getCHOGFileNames(querydir)
    total_in_range=0
    for querychog in queries:
        lat, lon = info.getQueryCHOGCoord(querychog)
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
  'checks': 1024,
  'trees': 1,
  'distance_type': 'kl',
  'descriptor': 'chog',
  'vote_method': 'highest',
  'confstring': '',
})
dbdump = os.path.join(maindir, "Research/collected_images/earthmine-new,culled/37.871955,-122.270829")
if __name__ == "__main__":
    querydir = os.path.join(maindir, 'Research/collected_images/query/%s/' % QUERY)
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
