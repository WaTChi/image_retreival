import os
import os.path
import shutil
import subprocess
import sys
import time

import groundtruthB
import groundtruthG
import groundtruthO
import groundtruthR
import groundtruthY
import info
import util

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
    queryimgpath = os.path.join(querydir, query + '.JPG')
    queryoutpath = os.path.join(resultsdir, query + ';query;gt' + str(match)  + ';' + dup + ';' + matchedimg + ';' + str(score) + '.jpg')
    shutil.copy(queryimgpath, queryoutpath)
    for matchedimg, score in ranked_matches:
        if score != topentry[1]:
            break
        matchimgpath = os.path.join(dbdump, '%s.jpg' % matchedimg)
        matchoutpath = os.path.join(resultsdir, query + ';match;gt' + str(match)  + ';' + dup + ';' + matchedimg + ';' + str(score) + '.jpg')
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

def query(querydir, querysift, dbdir, mainOutputDir, nClosestCells, copytopmatch, copy_top_n_percell=0):
    lat, lon = info.getQuerySIFTCoord(querysift)
    closest_cells = util.getclosestcells(lat, lon, dbdir)
    outputFilePaths = []
    cells_in_range = [(cell, dist) for cell, dist in closest_cells[0:nClosestCells] if dist < cellradius + ambiguity]
    if verbosity > 0:
        print "checking query: {0} \t against {1} \t cells".format(querysift, len(cells_in_range))
    for cell, dist in cells_in_range:
        if verbosity > 1:
            print "querying cell: {0}, distance: {1} with:{2}".format(cell, dist, querysift)
        outputFilePath = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(dist)  + ".res")
        outputFilePaths.append(outputFilePath)
        if  not os.path.exists(outputFilePath):
            subprocess.call(["/home/zhangz/workspace/Palantir/Debug/Palantir.thresh70,sparam1024,kd4", dbdir, cell, querydir, querysift, outputFilePath])
        if copy_top_n_percell > 0:
            outputDir = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(dist))
            copy_topn_results(os.path.join(dbdir, cell), outputDir, outputFilePath, 4)
#    combined = combine_until_dup(outputFilePaths, 1000)
    combined = combine_topn_votes(outputFilePaths, float('inf'))
    combined = weigh_cells(combined, dbdir, nClosestCells)
#    write_scores(querysift, combined, "/media/data/combined")
    [g, y, r, b, o] = check_topn_img(querysift, combined, topnresults)
    if copytopmatch:
        match = g or y or r or b or o
        copy_top_match(querydir, querysift.split('sift.txt')[0], combined, match)
    return [g, y, r, b, o]

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
        else:
            print "ERROR! weigh_cells has more overlap than it's supposed to: {0}, {1}".format(numoverlap, matchedimg)
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

def combine_until_dup(outputFilePaths, topn):
    #combines up to topn results from each query in outputFilePaths
    for n in range(1, topn + 1):
        combined = combine_topn_dup(outputFilePaths, n)
        if combined[0][1] > 1:
            break
    return combined

def combine_topn_dup(outputFilePaths, topn):
    #combines topn results from each query in outputFilePaths
    dupCount = {}
    for outputFilePath in outputFilePaths:
        for score, img in get_top_results(outputFilePath, topn):
            dupCount[img] = dupCount.get(img, 0) + 1
    dupCountLst = dupCount.items()
    dupCountLst.sort(key=lambda x: x[1], reverse=True)
    return dupCountLst

def combine_topn_votes(outputFilePaths, topn):
    #returns true if query if in topn of results
    dupCount = {}
    for outputFilePath in outputFilePaths:
        for score, img in get_top_results(outputFilePath, topn):
            dupCount[img] = dupCount.get(img, 0) + int(score)
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
        g += check_truth(querysift.split('sift')[0], entry[0], groundtruthG.matches)
        y += check_truth(querysift.split('sift')[0], entry[0], groundtruthY.matches)
        r += check_truth(querysift.split('sift')[0], entry[0], groundtruthR.matches)
        b += check_truth(querysift.split('sift')[0], entry[0], groundtruthB.matches)
        o += check_truth(querysift.split('sift')[0], entry[0], groundtruthO.matches)
    return [g > 0, y > 0, r > 0, b > 0, o > 0]

def check_top_score(querysift, dupCountLst):
    topCount = dupCountLst[0][1]
    g = 0
    y = 0
    r = 0
    b = 0
    o = 0
    count = 0
    for entry in dupCountLst:
        if entry[1] != topCount:
            print "# of imgs w/ top vote ({0}): {1}".format(topCount, count)
            break
        count += 1
        g += check_truth(querysift.split('sift')[0], entry[0], groundtruthG.matches)
        y += check_truth(querysift.split('sift')[0], entry[0], groundtruthY.matches)
        r += check_truth(querysift.split('sift')[0], entry[0], groundtruthR.matches)
        b += check_truth(querysift.split('sift')[0], entry[0], groundtruthB.matches)
        o += check_truth(querysift.split('sift')[0], entry[0], groundtruthO.matches)
    return [g > 0, y > 0, r > 0, b > 0, o > 0]

def characterize(querydir, dbdir, mainOutputDir, n, copytopmatch):
    start = time.time()
    if not os.path.exists(mainOutputDir):
        os.makedirs(mainOutputDir)
    if copytopmatch:
        if os.path.exists(resultsdir):
            shutil.rmtree(resultsdir)
        os.makedirs(resultsdir)
    files = util.getSiftFileNames(querydir)
    g_count = 0
    y_count = 0
    r_count = 0
    b_count = 0
    o_count = 0
    count = 0
    for queryfile in files:
        count += 1
        [g, y, r, b, o] = query(querydir, queryfile, dbdir, mainOutputDir, n, copytopmatch)
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
    print "g:{0} y:{1} r:{2} b:{3} o:{4} = {5}, out of {6}={7}%".format(g_count, y_count, r_count, b_count, o_count, total_count, count, match_rate)

#def convertGroundTruth(gt,file):
#    for key in gt.keys:


cellradius = 236.6
ambiguity = 50
ncells = 7 #if ambiguity<100, 7 is max possible by geometry
topnresults = 1
verbosity = 0
copytopmatch = False
resultsdir = '/media/data/topmatches'
maindir = "/home/zhangz/.gvfs/data on 128.32.43.40"
dbdump = os.path.join(maindir, "Research/collected_images/earthmine-new,culled/37.871955,-122.270829")
if __name__ == "__main__":
    querydir = os.path.join(maindir, 'query3/')
    dbdir = os.path.join(maindir, 'Research/cellsg=100,r=d=236.6/')
#    matchdir = os.path.join(maindir, 'Research/results(query3)/matchescells(g=100,r=d=236.6),query3,kdtree4,threshold=70k,searchparam=2048')
    matchdir = os.path.join(maindir, 'Research/results(query3)/matchescells(g=100,r=d=236.6),query3,kdtree4,threshold=70k,searchparam=1024')
    if len(sys.argv) > 4:
        print "USAGE: {0} QUERYDIR DBDIR OUTPUTDIR".format(sys.argv[0])
        sys.exit()
    elif len(sys.argv) == 4:
        querydir = sys.argv[1]
        dbdir = sys.argv[2]
        matchdir = sys.argv[3]
    for n in [50]:
        ambiguity = n
        print "ambiguity:{0}".format(n)
        for n2 in [2]:
            topnresults = n2
            print "\t top {0} results".format(n2)
            characterize(querydir, dbdir, matchdir, ncells, copytopmatch)

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