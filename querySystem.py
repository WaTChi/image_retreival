import os
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

def check_top(query_str, groundTruth_dict, outputFilePath, n=1):
    if  not os.path.exists(outputFilePath):
        raise OSError("{p} does not exist.".format(p=outputFilePath))
    file = open(outputFilePath)
    top_results = []
    line = file.next()
    score1, img = parse_result_line(line)
    top_results.append(img)
    uniquevotecount = 0
    for line in file:
        score2, img = parse_result_line(line)
        if score2 != score1:
            score1 = score2
            uniquevotecount += 1
            if uniquevotecount >= n:
                break
        top_results.append(img)
    for img in top_results:
        if check_truth(query_str.split('sift')[0], img, groundTruth_dict):
            return True
    return False

def copy_query(querydir, querysift, outdir):
    if  not os.path.exists(outdir):
        os.makedirs(outdir)
    img = '%s.JPG' % querysift.split('sift')[0]
    shutil.copy(os.path.join(querydir, img), outdir)

def copy_results(imgdir, outdir, filepath):
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

def query(querydir, querysift, dbdir, mainOutputDir, nClosestCells):
    lat, lon = info.getQuerySIFTCoord(querysift)
    closest_cells = util.getclosestcells(lat, lon, dbdir)
    outputFilePaths = []
    for cell, dist in closest_cells[0:nClosestCells]:
        if dist > 286.6:
            continue
        print "querying cell: {0}, distance: {1} with:{2}".format(cell, dist, querysift)
        outputFilePath = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(dist)  + ".res")
        #outputFilePath = os.path.join(mainOutputDir, querysift + ',' + cell + ".res")
        outputFilePaths.append(outputFilePath)
        if  not os.path.exists(outputFilePath):
            subprocess.call(["/home/zhangz/workspace/Palantir/Debug/Palantir.thresh70,sparam1024,kd4", dbdir, cell, querydir, querysift, outputFilePath])
#        outputDir = os.path.join(mainOutputDir, querysift + ',' + cell + ',' + str(dist))
#        copy_results(os.path.join(dbdir, cell), outputDir, outputFilePath)
#    write_topn_img(querysift, outputFilePaths, "/media/data/query3top2", 2)
    return check_topn_img(querysift, outputFilePaths, 1)

def write_topn_img(querysift, outputFilePaths, outdir, topnres=10, topn=float("inf")):
    if  not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = open(os.path.join(outdir, querysift + ".top"), 'w')
    dupCountLst = combine_topn_votes(outputFilePaths, topn)
    for entry in dupCountLst[0:topnres]:
        outfile.write(entry[0])
        outfile.write('\n')

def get_top_results(outputFilePath, n):
    file = open(outputFilePath)
    top_results = []
    i = 0
    for line in file:
        if i > n:
            break
        score, img = parse_result_line(line)
        top_results.append((score, img))
        i += 1
    file.close()
    return top_results

def combine_topn_votes(outputFilePaths, topn=float("inf")):
    #combines topn results from each query in outputFilePaths
    #returns true if query if in topnres
    dupCount = {}
    for outputFilePath in outputFilePaths:
        for score, img in get_top_results(outputFilePath, topn):
            dupCount[img] = dupCount.get(img, 0) + int(score)
    dupCountLst = dupCount.items()
    dupCountLst.sort(key=lambda x: x[1])
    dupCountLst.reverse()
    return dupCountLst

def check_topn_img(querysift, outputFilePaths, topnres=1, topn=float("inf")):
    dupCountLst = combine_topn_votes(outputFilePaths, topn)
    g = 0
    y = 0
    r = 0
    b = 0
    o = 0
    for entry in dupCountLst[0:topnres]:
        # print "votes: {0}".format(entry[1])
            g += check_truth(querysift.split('sift')[0], entry[0], groundtruthG.matches)
            y += check_truth(querysift.split('sift')[0], entry[0], groundtruthY.matches)
            r += check_truth(querysift.split('sift')[0], entry[0], groundtruthR.matches)
            b += check_truth(querysift.split('sift')[0], entry[0], groundtruthB.matches)
            o += check_truth(querysift.split('sift')[0], entry[0], groundtruthO.matches)
    return [g > 0, y > 0, r > 0, b > 0, o > 0]

    
def check_top_combined_vote(querysift, outputFilePaths, combinetopn=float("inf")):
    dupCountLst = combine_topn_votes(outputFilePaths, combinetopn)
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
    
def check_top_dup(querysift, outputFilePaths, combinetopn=5):
    dupCount = {}
    for outputFilePath in outputFilePaths:
        for score, img in get_top_results(outputFilePath, combinetopn):
            dupCount[img] = dupCount.get(img, 0) + 1
    dupCountLst = dupCount.items()
    dupCountLst.sort(key=lambda x: x[1])
    dupCountLst.reverse()
    topCount = dupCountLst[0][1]
    g = 0
    y = 0
    r = 0
    b = 0
    o = 0
    count = 0
    for entry in dupCountLst:
        if entry[1] != topCount:
            print "# of imgs w/ top dup ({0}): {1}".format(topCount, count)
            break
        count += 1
        g += check_truth(querysift.split('sift')[0], entry[0], groundtruthG.matches)
        y += check_truth(querysift.split('sift')[0], entry[0], groundtruthY.matches)
        r += check_truth(querysift.split('sift')[0], entry[0], groundtruthR.matches)
        b += check_truth(querysift.split('sift')[0], entry[0], groundtruthB.matches)
        o += check_truth(querysift.split('sift')[0], entry[0], groundtruthO.matches)
    return [g > 0, y > 0, r > 0, b > 0, o > 0]

def characterize(querydir, dbdir, mainOutputDir, n=1):
    start = time.time()
    if not os.path.exists(mainOutputDir):
        try:
            os.makedirs(mainOutputDir)
        except Exception:
            print "Error making directory {0} ...quitting...".format(mainOutputDir)
            return
    baddir = "/media/data/Research/bad"
    gooddir = "/media/data/Research/good"
    if  os.path.exists(baddir):
        shutil.rmtree(baddir)
    if  os.path.exists(gooddir):
        shutil.rmtree(gooddir)
    os.makedirs(baddir)
    os.makedirs(gooddir)
    files = util.getSiftFileNames(querydir)
    g_count = 0
    y_count = 0
    r_count = 0
    b_count = 0
    o_count = 0
    count = 0
    for queryfile in files:
        count += 1
        print "checking query: {0}".format(queryfile)
        [g, y, r, b, o] = query(querydir, queryfile, dbdir, mainOutputDir, n)
        if not (g or y or r or b or o):
            copy_query(querydir, queryfile, baddir)
#        else:
#            copy_query(querydir, queryfile, gooddir)
        if g:
            g_count += 1
            print "G match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        elif y:
            y_count += 1
            print "Y match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        elif r:
            r_count += 1
            print "R match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        elif b:
            b_count += 1
            print "B match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        elif o:
            o_count += 1
            print "O match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
        else:
            print "No match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
    end = time.time()
    elapsed = end - start
    print "total time:{0}, avg time:{1}".format(elapsed, elapsed / count)

    total_count = g_count + y_count + r_count + b_count + o_count
    match_rate = float(total_count) / count
    print "g:{0} y:{1} r:{2} b:{3} o:{4} = {5}, out of {6}={7}%".format(g_count, y_count, r_count, b_count, o_count, total_count, count, match_rate)




def print_cells(querydir, dbdir):
    for querysift in util.getSiftFileNames(querydir):
        lat, lon = info.getQuerySIFTCoord(querysift)
        print querysift
        util.printclosestcells(lat, lon, dbdir, 336.6)

if __name__ == "__main__":

    n = 1
#    querydir = "E:\Research\collected_images\query\query3"
#    dbdir = "E:\Research\cellsg=100,r=d=236.6"
#    matchdir="E:\cellsg=100,r=d=236.6,query3kdtree1matches,threshold=70k"

#    querydir = "/media/data/Research/query3/"
#    dbdir = "/media/data/Research/cellsg=100,r=d=236.6/"
    matchdir = "/media/data/Research/matchescells(g=100,r=d=236.6),query3,kdtree1,threshold=70k,searchparam256/"

##    matchdir = "/home/zhangz/.gvfs/data on 128.32.43.40/Research/results(query3)/matchescells(g=100,r=d=236.6),query3,kdtree1,threshold=70k/"
    querydir = "/home/zhangz/.gvfs/data on 128.32.43.40/query3/"
    dbdir = "/home/zhangz/.gvfs/data on 128.32.43.40/Research/cellsg=100,r=d=236.6/"
#    matchdir = "/home/zhangz/.gvfs/data on 128.32.43.40/Research/results(query3)/matchescells(g=100,r=d=236.6),query3,kdtree4,threshold=70k,searchparam=1024"
    matchdir = "/home/zhangz/.gvfs/data on 128.32.43.40/Research/results(query3)/matchescells(g=100,r=d=236.6),query3,kdtree1,threshold=70k"

    if len(sys.argv) > 4:
        print "USAGE: {0} QUERYDIR DBDIR OUTPUTDIR".format(sys.argv[0])
        sys.exit()
    elif len(sys.argv) == 4:
        querydir = sys.argv[1]
        dbdir = sys.argv[2]
        matchdir = sys.argv[3]
    #util.writeQueryCoords(querydir, "E:\query3.txt")
    characterize(querydir, dbdir, matchdir, n)
