import os
import os.path
import shutil
import sys
import time

import groundtruthB
import groundtruthG
import groundtruthO
import groundtruthR
import groundtruthY
import info
import query2Groundtruth
import util
import query
import datetime

HOME = os.path.expanduser('~')

def parse_result_line(line):
    score = line.split('\t')[0]
    img = line.split('\t')[1].split('sift')[0]
    return score, img

def check_truth(query_str, result_str, groundTruth_dict):
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
        shutil.copy(os.path.join(imgdir, ('%s.jpg' % img)), os.path.join(outdir, '%s-%s.jpg' % (score, img)))

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

    queryimgpath = os.path.join(querydir, query + '.JPG')
    queryoutpath = os.path.join(resultsdir, query + ';query;gt' + str(match)  + ';' + dup + ';' + matchedimg + ';' + str(score) + ';' + str(distance) + '.jpg')
    shutil.copy(queryimgpath, queryoutpath)
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

def run_query(newlat, newlon, querydir, querysift, dbdir, mainOutputDir, nClosestCells, copytopmatch, closest_cells, params, copy_top_n_percell=0):
    outputFilePaths = []
    cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:nClosestCells] if dist < cellradius + ambiguity+matchdistance]
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
    query.run_parallel(dbdir, [c for c,d in cells_in_range], querydir, querysift, outputFilePaths, params)
#     end query
    if copy_top_n_percell > 0:
        for cell, dist in cells_in_range:
            latcell, loncell = cell.split(',')
            latcell = float(latcell)
            loncell = float(loncell)
            actualdist = info.distance(latquery, lonquery, latcell, loncell)
            outputFilePath = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(actualdist)  + ".res")
            outputDir = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(actualdist))
            copy_topn_results(os.path.join(dbdir, cell), outputDir, outputFilePath, 4)
##    combined = combine_until_dup(outputFilePaths, 1000)
#    combined = combine_topn_votes(outputFilePaths, float('inf'))
##    combined = filter_in_range(newlat, newlon, combined, querysift)
##    write_scores(querysift, combined, "/media/data/combined")
    results = {}
    for n in topnresults:
#        result = check_topn_img(querysift, combined, n)
#        results[n] = reduce(lambda x,y: x or y, result)
        results[n]=False
#    if copytopmatch:
#        copy_top_match(querydir, querysift.split('sift.txt')[0], combined, results[1])
    return results

def filter_in_range(qlat, qlon, ranked_matches, querysift):
    #qlat, qlon = info.getQuerySIFTCoord(querysift)
    filtered_matches = []
    for matchedimg, score in ranked_matches:
        clat = float(matchedimg.split(",")[0])
        clon = float(matchedimg.split(",")[1][0:-5])
        distance = info.distance(qlat, qlon, clat, clon)
        if distance < ambiguity+matchdistance:
            filtered_matches.append((matchedimg, score))
#    filtered_matches.sort(key=lambda x: x[1], reverse=True)
    return filtered_matches


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
    dupCountLst.sort(key=lambda x: x[1], reverse=True)
    return dupCountLst

def check_topn_img(querysift, dupCountLst, topnres=1):
    g = 0
    y = 0
    r = 0
    b = 0
    o = 0
    for entry in dupCountLst[0:topnres]:
#        g += check_truth(querysift.split('sift')[0], entry[0], query2Groundtruth.matches)

        g += check_truth(querysift.split('sift')[0], entry[0], groundtruthG.matches)
        y += check_truth(querysift.split('sift')[0], entry[0], groundtruthY.matches)
        r += check_truth(querysift.split('sift')[0], entry[0], groundtruthR.matches)
        b += check_truth(querysift.split('sift')[0], entry[0], groundtruthB.matches)
        o += check_truth(querysift.split('sift')[0], entry[0], groundtruthO.matches)
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

def characterize_fuzzy(querydir, dbdir, mainOutputDir, ncells, copytopmatch, params):
    start = time.time()
    if not os.path.exists(mainOutputDir):
        os.makedirs(mainOutputDir)
    if copytopmatch:
        if os.path.exists(resultsdir):
            shutil.rmtree(resultsdir)
        os.makedirs(resultsdir)
    files = util.getSiftFileNames(querydir)
    results={}
    for n in topnresults:
        results[n]=0
    count = 0
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

#def get_num_imgs_in_range(range, celldir, querydir=HOME + '/.gvfs/data on 128.32.43.40/query3/'):
#    queries = util.getSiftFileNames(querydir)
#    total_in_range=0
#    for querysift in queries:
#        lat, lon = info.getQuerySIFTCoord(querysift)
#        cell, dist = util.getclosestcell(lat, lon, celldir)
#        numinrange = util.getNumJPGInRange(lat, lon, os.path.join(celldir,cell), range)
#        total_in_range+=numinrange
#    print "average # of imgs within {0} meters of query: {1}".format(range, float(total_in_range)/len(queries))

cellradius = 236.6
ambiguity = 0
matchdistance = 25
ncells =  7  #if ambiguity+matchdistance<100, 7 is max possible by geometry
topnresults = [1, 2,5,10]
verbosity = 0
copytopmatch = False
resultsdir = '/media/data/topmatches'
maindir = HOME + "/shiraz"
fuzzydir = os.path.join(maindir, 'query3_fuzzy2/')
dbdump = os.path.join(maindir, "Research/collected_images/earthmine-new,culled/37.871955,-122.270829")
params = query.PARAMS_DEFAULT.copy()
params.update({
  'checks': 1024,
  'algorithm': 'kdtree',
  'trees': 4,
  'vote_method': 'highest',
  'confstring': '',
})
if __name__ == "__main__":
    querydir = os.path.join(maindir, 'query3/')
    dbdir = os.path.join(maindir, 'Research/cellsg=100,r=d=236.6/')
    #for ahvaz
#    matchdir = HOME+'/.gvfs/sftp for zhangz on 128.32.43.34/home/zhangz/results(%s)/matchescells(g=100,r=d=236.6),%s,%s' % ('query3', 'query3', query.searchtype(params))
    #for gorgan
    matchdir = HOME+'/results(%s)/matchescells(g=100,r=d=236.6),%s,%s' % ('query3', 'query3', query.searchtype(params))
    #for q2
#    matchdir = os.path.join(maindir, 'Research/results(%s)/matchescells(g=100,r=d=236.6),%s,%s' % ('query2', 'query2', query.searchtype(params)))
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
        #get_num_imgs_in_range(ambiguity+matchdistance, dbdir)
        characterize_fuzzy(querydir, dbdir, matchdir, ncells, copytopmatch, params)
        print "ambiguity:{0}".format(n)
        print querydir
        print dbdir
        print matchdir


#def combine_until_dup(outputFilePaths, topn):
#    #combines up to topn results from each query in outputFilePaths
#    for n in range(1, topn + 1):
#        combined = combine_topn_dup(outputFilePaths, n)
#        if combined[0][1] > 1:
#            break
#    return combined
#
#def combine_topn_dup(outputFilePaths, topn):
#    #combines topn results from each query in outputFilePaths
#    dupCount = {}
#    for outputFilePath in outputFilePaths:
#        for score, img in get_top_results(outputFilePath, topn):
#            dupCount[img] = dupCount.get(img, 0) + 1
#    dupCountLst = dupCount.items()
#    dupCountLst.sort(key=lambda x: x[1], reverse=True)
#    return dupCountLst
#
#def check_top_score(querysift, dupCountLst):
#    topCount = dupCountLst[0][1]
#    g = 0
#    y = 0
#    r = 0
#    b = 0
#    o = 0
#    count = 0
#    for entry in dupCountLst:
#        if entry[1] != topCount:
#            print "# of imgs w/ top vote ({0}): {1}".format(topCount, count)
#            break
#        count += 1
#        g += check_truth(querysift.split('sift')[0], entry[0], groundtruthG.matches)
#        y += check_truth(querysift.split('sift')[0], entry[0], groundtruthY.matches)
#        r += check_truth(querysift.split('sift')[0], entry[0], groundtruthR.matches)
#        b += check_truth(querysift.split('sift')[0], entry[0], groundtruthB.matches)
#        o += check_truth(querysift.split('sift')[0], entry[0], groundtruthO.matches)
#    return [g > 0, y > 0, r > 0, b > 0, o > 0]

#def copy_query(querydir, querysift, outdir):
#    if  not os.path.exists(outdir):
#        os.makedirs(outdir)
#    img = '%s.JPG' % querysift.split('sift')[0]
#    shutil.copy(os.path.join(querydir, img), outdir)

#def check_top(query_str, groundTruth_dict, outputFilePath, n=1):
#    if  not os.path.exists(outputFilePath):
#        raise OSError("{p} does not exist.".format(p=outputFilePath))
#    file = open(outputFilePath)
#    top_results = []
#    line = file.next()
#    score1, img = parse_result_line(line)
#    top_results.append(img)
#    uniquevotecount = 0
#    for line in file:
#        score2, img = parse_result_line(line)
#        if score2 != score1:
#            score1 = score2
#            uniquevotecount += 1
#            if uniquevotecount >= n:
#                break
#        top_results.append(img)
#    for img in top_results:
#        if check_truth(query_str.split('sift')[0], img, groundTruth_dict):
#            return True
#    return False

#def weigh_cells(ranked_matches, dbdir, nClosestCells):
#    weighted_matches = []
#    for matchedimg, score in ranked_matches:
#        lat = float(matchedimg.split(",")[0])
#        lon = float(matchedimg.split(",")[1][0:-5])
#        closest_cells = util.getclosestcells(lat, lon, dbdir)
#        cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:nClosestCells] if dist < cellradius]
#        numoverlap = len(cells_in_range)
#        if numoverlap == 3:
#            weighted_matches.append((matchedimg, score * 4))
#        elif numoverlap == 4:
#            weighted_matches.append((matchedimg, score * 3))
##        else:
##            print "ERROR! weigh_cells has more overlap than it's supposed to: {0}, {1}".format(numoverlap, matchedimg)
#    weighted_matches.sort(key=lambda x: x[1], reverse=True)
#    return weighted_matches