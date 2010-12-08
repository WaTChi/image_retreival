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
    img = line.split('\t')[1].split('.jpg')[0]
    return score, img

def check_truth(query_str, result_str, groundTruth_dict):
    return result_str in groundTruth_dict[query_str]

featureExtrator = "/home/zhangz/Downloads/libpmk-2.5/libpmk_features/tools/images-to-pointsets.out"
pyramid_builder = "/home/zhangz/Downloads/libpmk-2.5/libpmk2/tools/clusters-to-pyramids.out"

def query(querydir, querysift, dbdir, mainOutputDir, nClosestCells):
    lat, lon = info.getQuerySIFTCoord(querysift)
    closest_cells = util.getclosestcells(lat, lon, dbdir)
    outputFilePaths = []

    queryJPGpath = os.path.join(querydir,querysift.split('sift.txt')[0]+'.JPG')
    queryPSLpath =os.path.join(querydir, querysift.split('sift.txt')[0]+'.psl')
    if not os.path.exists(queryPSLpath):
        print "making psl"
        subprocess.call([featureExtrator, queryPSLpath, queryJPGpath], cwd='/home/zhangz/Downloads/libpmk-2.5/libpmk_features/tools')

    for cell, dist in closest_cells[0:nClosestCells]:
        if dist > 150:
            continue

        querymrh =os.path.join(querydir, querysift.split('sift.txt')[0]+','+cell+'.mrh')
        if not os.path.exists(querymrh):
            print "making mrh"
            hcpath = os.path.join(dbdir, cell+'.hc')
            subprocess.call([pyramid_builder, queryPSLpath, hcpath, querymrh], cwd='/home/zhangz/Downloads/libpmk-2.5/libpmk2/tools/')

        print "querying cell: {0}, distance: {1} with:{2}".format(cell, dist, querysift)
        dbpath = os.path.join(dbdir, cell+'.mrh')
        
        outputFilePath = os.path.join(mainOutputDir, querysift.split('sift.txt')[0] + ',' + cell + ',' + str(dist)  + ".pyres")
        outputFilePaths.append(outputFilePath)
        if  not os.path.exists(outputFilePath):
            subprocess.call(["/home/zhangz/Downloads/libpmk-2.5/libpmk2/tools/horus.out", querymrh, dbpath, outputFilePath])
        outputDir = os.path.join(mainOutputDir, querysift.split('sift.txt')[0] + ',' + cell + ',' + str(dist))
        print outputFilePath
        copy_results(os.path.join(dbdir, cell), outputDir, outputFilePath)
    if len(outputFilePaths)==0:
        return [0,0,0,0,0]
    return check_top_dup(querysift, outputFilePaths)



def copy_results(imgdir, outdir, filepath):
    if  os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    results = get_top_results(filepath)
    for score, img in results[0:4]:
        shutil.copy(os.path.join(imgdir, ('%s.jpg' % img)), os.path.join(outdir, '%s-%s.jpg' % (score, img)))

    for score, img in results[-4:len(results)]:
        shutil.copy(os.path.join(imgdir, ('%s.jpg' % img)), os.path.join(outdir, '%s-%s.jpg' % (score, img)))

def getfilemapping(mappingfilepath):
    mapping = []
    file = open(mappingfilepath)
    for line in file:
        img = line.split('\t')[1].split('.jpg')[0]
        mapping.append(img)
    return mapping

def get_top_results(outputFilePath, n=float('inf')):
    file = open(outputFilePath)

    lat = outputFilePath.split("DSC_")[1].split(',')[3]
    lon = outputFilePath.split("DSC_")[1].split(',')[4]

    mapping = getfilemapping("/media/data/Research/cellsg=50,r=100,d=86.6/" + lat+','+lon+'.map')

    results = []
    for line in file:
        score, idx = parse_result_line(line)
        img = mapping[int(idx)]
        results.append((float(score), img))
    file.close()
    results.sort(key=lambda x: x[0])
    results.reverse()
    if n==float("inf"):
        return results
    return results[0:n]

def combine_topn_votes(outputFilePaths, topn=float("inf")):
    #combines topn results from each query in outputFilePaths
    #returns true if query if in topnres
    dupCount = {}
    for outputFilePath in outputFilePaths:
        for score, img in get_top_results(outputFilePath, topn):
            dupCount[img] = dupCount.get(img, 0) + float(score)
    dupCountLst = dupCount.items()
    dupCountLst.sort(key=lambda x: x[1])
    #dupCountLst.reverse()
    return dupCountLst

def check_topn_img(querysift, outputFilePaths, topnres=1, topn=float("inf")):
    dupCountLst = combine_topn_votes(outputFilePaths, topn)
    g = 0
    y = 0
    r = 0
    b = 0
    o = 0
    print "first match: {0}".format(dupCountLst[0])
    print "last match: {0}".format(dupCountLst[-1])
    for entry in dupCountLst[0:topnres]:
        # print "votes: {0}".format(entry[1])
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





if __name__ == "__main__":

    n = 7
    querydir = "/media/data/Research/query3/"
    dbdir = "/media/data/Research/cellsg=50,r=100,d=86.6/"
    matchdir = "/media/data/Research/matchescells(g=50,r=100,d=86.6),query3,pmk/"

    if len(sys.argv) > 4:
        print "USAGE: {0} QUERYDIR DBDIR OUTPUTDIR".format(sys.argv[0])
        sys.exit()
    elif len(sys.argv) == 4:
        querydir = sys.argv[1]
        dbdir = sys.argv[2]
        matchdir = sys.argv[3]
    characterize(querydir, dbdir, matchdir, n)
