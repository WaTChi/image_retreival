# The query system.
# Collection of functions that act on contexts.

import os
import os.path
import query
from config import INFO
import posit
import pnp

try:
    if 'NUM_THREADS' in os.environ:
        NUM_THREADS = int(os.environ['NUM_THREADS'])
    else:
        import multiprocessing
        NUM_THREADS = multiprocessing.cpu_count()
except:
    import multiprocessing
    NUM_THREADS = multiprocessing.cpu_count()

import time

from config import *
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
import util_cs188

from multiprocessing import Pool, cpu_count

class MultiprocessExecution:
  pool = None

  def __enter__(self):
    MultiprocessExecution.pool = Pool(cpu_count())

  def __exit__(self, *args):
    print "Waiting for background jobs to finish..."
    MultiprocessExecution.pool.close()
    MultiprocessExecution.pool.join()
    MultiprocessExecution.pool = None
    print "All processes done."

class LocationOutOfRangeError(Exception):
    """Raised when there are no cells near query location"""
    pass

def check_truth(query_str, result_str, groundTruth_dict):
    return result_str in groundTruth_dict[query_str]

def skew_location(C, Q):
    length = 2*C.ambiguity
    points = []
    corner = info.moveLocation(Q.sensor_lat, Q.sensor_lon, (2**.5)*C.ambiguity, -45)
    for i in range(length+1):
        row = info.moveLocation(corner[0], corner[1], i, 180)
        for j in range(length+1):
            point = info.moveLocation(row[0], row[1], j, 90)
            if info.distance(Q.sensor_lat, Q.sensor_lon, point[0], point[1]) <= C.ambiguity:
                points.append(point)
    return points

def load_location(C, Q):
    file = open(os.path.join(C.fuzzydir, Q.name)+'.fuz')
    points = []
    for line in file:
        lat, lon = line.strip().split('\t')
        lat = float(lat)
        lon = float(lon)
        points.append((lat,lon))
    return points

def derive_key(closest_cells, name):
    return (name,) + tuple(sorted(map(lambda (cell, dist): cell, closest_cells)))

def distance_sort(C, Q, matches):
  def extract_key(match):
    line = match[0]
    lat, lon = getlatlonfromdbimagename(C, line)
    return (match[1], -info.distance(lat, lon, Q.query_lon, Q.query_lat))
  return sorted(matches, key=extract_key, reverse=True)

cache = {}
def match(C, Q):
    # compute closest cells
    closest_cells = util.getclosestcells(Q.query_lat, Q.query_lon, C.dbdir)
    cells_in_range = [(cell, dist)
      for cell, dist in closest_cells[0:C.ncells]
        if dist < C.cellradius + C.ambiguity + C.matchdistance]

    if not cells_in_range:
        raise LocationOutOfRangeError

    # cache for fuzz runs
    if C.cacheEnable:
        key = derive_key(cells_in_range, Q.siftname)
        if key in cache:
            return cache[key]

    # compute output file paths for the cells
    outputFilePaths = []
    for cell, dist in cells_in_range:
        latcell, loncell = cell.split(',')
        latcell = float(latcell)
        loncell = float(loncell)
        actualdist = info.distance(Q.query_lat, Q.query_lon, latcell, loncell)
        outputFilePath = os.path.join(C.matchdir, Q.siftname + ',' + cell + ',' + str(actualdist)  + ".res")
        outputFilePaths.append(outputFilePath)

    # start query
    query.run_parallel(C, Q, [c for c,d in cells_in_range], outputFilePaths, NUM_THREADS)

    # combine results
    comb_matches = corr.combine_matches(outputFilePaths)
    ranked = distance_sort(C, Q, \
      combine_ransac(comb_matches, C))
#    print cells_in_range
#    for i,r in enumerate(ranked):
#      print i, r

    # top 1
    stats = check_topn_img(C, Q, ranked, 1)

    # compute homography and draw images maybe
    if C.compute_hom:
      if MultiprocessExecution.pool:
        MultiprocessExecution.pool.apply_async(compute_hom, [C.pickleable(), Q, ranked, comb_matches])
      else:
        compute_hom(C, Q, ranked, comb_matches)

    # return statistics and top result
    matchedimg = ranked[0][0]
    matches = comb_matches[matchedimg + 'sift.txt']
    if C.cacheEnable:
        cache[key] = (stats, matchedimg, matches, ranked)
    if C.match_callback:
        C.match_callback(C, Q, stats, matchedimg, ranked, cells_in_range)

    # done
    return stats, matchedimg, matches, ranked

def getlatlonfromdbimagename(C, dbimg):
    if C.QUERY == 'emeryville':
        return 0,0
    clat = float(dbimg.split(",")[0])
    clon = float(dbimg.split(",")[1][0:-5])
    return clat, clon
    
# top entry is ranked_matches[0], etc
def compute_hom(C, Q, ranked_matches, comb_matches):
    if C.put_into_dirs:
        udir = os.path.join(C.resultsdir, Q.name)
    else:
        udir = C.resultsdir
    if not os.path.exists(udir):
        os.makedirs(udir)
    i = 0
    data = {}
    for matchedimg, score in ranked_matches[:C.max_matches_to_analyze]:
        i += 1
        if C.stop_on_homTrue and data.get('success'):
            break # we are done (found good homography)
        clat, clon = getlatlonfromdbimagename(C, matchedimg)
        matchimgpath = os.path.join(C.dbdump, '%s.jpg' % matchedimg)
        if not os.path.exists(matchimgpath):
            matchimgpath = os.path.join(C.dbdump, '%s.JPG' % matchedimg)
            assert os.path.exists(matchimgpath)
        match = any(check_img(C, Q, ranked_matches[i-1]))

        # rematch for precise fit
        db_matches = comb_matches[matchedimg + 'sift.txt']
        matches = db_matches
        matchsiftpath = os.path.join(C.dbdump, matchedimg + 'sift.txt')
        matches = corr.rematch(C, Q, matchsiftpath)

        # concat db matches
        matches.extend(db_matches)

        # find homography
        rsc_matches, H, inliers = corr.find_corr(matches, hom=True, ransac_pass=True, data=data)
        rsc_inliers = np.compress(inliers, rsc_matches).tolist()
        u = corr.count_unique_matches(rsc_inliers)

#        # compute pose [experimental]
#        if data.get('success'):
#          corr.compute_pose(C, rsc_matches, matchimgpath, matchsiftpath)

        if C.drawtopcorr:
          # draw picture
          matchoutpath = os.path.join(udir, Q.name + ';match' + str(i) + ';gt' + str(match)  + ';hom' + str(data.get('success')) + ';uniq=' + str(u) + ';inliers=' + str(float(sum(inliers))/len(matches)) + ';' + matchedimg + '.jpg')
          try:
            corr.draw_matches(C, Q, matches, rsc_matches, H, inliers, matchimgpath, matchoutpath)
          except IOError, e:
            INFO(e)

        if C.dump_hom:
          if C.put_into_dirs:
              identifier = str(i)
          else:
              identifier = Q.name + ':' + str(i)
          H = np.matrix(np.asarray(H))
          with open(os.path.join(udir, 'homography%s.txt' % identifier), 'w') as f:
              print >> f, H
          np.save(os.path.join(udir, 'inliers%s.npy' % identifier), rsc_inliers)

        ### POSIT ###
        if C.do_posit:
            posit.do_posit(C, Q, rsc_inliers, matchsiftpath, matchimgpath)

        ### Perspective N-Point Problem ###
        if C.solve_pnp:
            pnp.solve(C, Q, rsc_inliers, matchsiftpath, matchimgpath)
        

def combine_ransac(counts, C):
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    filtered = {}
    bound = -1
    num_filt = 0
    for siftfile, matches in sorted_counts:
      siftfile = siftfile[:-8]
      if num_filt > C.ransac_min_filt and (len(matches) < bound or num_filt > C.ransac_max_filt):
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

def check_img(C, Q, entry):
    g,y,r,b,o = 0,0,0,0,0
    if C.QUERY == 'query1':
        g += check_truth(Q.name, entry[0], query1GroundTruth.matches)
    elif C.QUERY == 'query3':
        g += check_truth(Q.name, entry[0], groundtruthG.matches)
        y += check_truth(Q.name, entry[0], groundtruthY.matches)
        r += check_truth(Q.name, entry[0], groundtruthR.matches)
        b += check_truth(Q.name, entry[0], groundtruthB.matches)
        o += check_truth(Q.name, entry[0], groundtruthO.matches)
    elif C.QUERY == 'query2':
        g += check_truth(Q.name, entry[0], query2Groundtruth.matches)
    elif C.QUERY == 'query4':
        g += check_truth(Q.name, entry[0], query4GroundTruth.matches)
    else:
        return [0,0,0,0,0]
    return [g > 0, y > 0, r > 0, b > 0, o > 0]

def check_topn_img(C, Q, dupCountLst, topnres=1):
    record = [0]*5
    for entry in dupCountLst[0:topnres]:
        new = check_img(C, Q, entry)
        record = map(lambda a,b: a + b, record, new)
    return map(bool, record)

def dump_combined_matches(C, Q, stats, matchedimg, matches, cells_in_range):
    # For Aaron's analysis
    table = {}
    for line in open(os.path.join(C.dbdir, 'cellmap.txt')):
        a, b = line.split()
        table[b] = int(a)
    def cellsetstr(cells):
        cells = sorted(map(lambda (cell, dist): str(table[cell]), cells))
        return '-'.join(cells)
    outputFilePath = os.path.join(C.matchdir, 'fuzz2', Q.siftname + ',combined,' + cellsetstr(cells_in_range) + ".res")
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

def characterize(C):
    INFO("matchdir=%s" % C.matchdir)
    start = time.time()
    results = util_cs188.Counter()
    C.initdirs()
    count = 0
    g_count = 0
    y_count = 0
    r_count = 0
    b_count = 0
    o_count = 0
    for Q in C.iter_queries():
        print '-- query', Q.name, '--'
        for loc in C.locator_function(C, Q):
            Q.setQueryCoord(*loc)
            count += 1
            try:
                [g, y, r, b, o], matchedimg, matches, combined = match(C, Q)
            except LocationOutOfRangeError:
                INFO('Exception: location out of cell range')
                continue
            # compile statistics
            # top n
            for n in C.topnresults:
                result = check_topn_img(C, Q, combined, n)
                results[n] += reduce(lambda x,y: x or y, result)
            if count % C.print_per == 0 and C.verbosity > 0:
                INFO('speed is %f' % ((time.time()-start)/count))
                for n in C.topnresults:
                    print "matched {0}\t out of {1}\t in the top {2}\t amb: {3}, ncells:{4}".format(results[n], count, n, C.ambiguity, C.ncells)

            if g:
                g_count += 1
                if C.verbosity > 0 and count % C.print_per == 0:
                    print "G match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            elif y:
                y_count += 1
                if C.verbosity > 0 and count % C.print_per == 0:
                    print "Y match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            elif r:
                r_count += 1
                if C.verbosity > 0 and count % C.print_per == 0:
                    print "R match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            elif b:
                b_count += 1
                if C.verbosity > 0 and count % C.print_per == 0:
                    print "B match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            elif o:
                o_count += 1
                if C.verbosity > 0 and count % C.print_per == 0:
                    print "O match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)
            else:
                if C.verbosity > 0 and count % C.print_per == 0:
                    print "No match-g:{0} y:{1} r:{2} b:{3} o:{4} out of {5}".format(g_count, y_count, r_count, b_count, o_count, count)

    end = time.time()
    elapsed = end - start
    if not count:
      print "Query set empty!"
      return results, 0
    if C.verbosity > 0:
        print "total time:{0}, avg time:{1}".format(elapsed, elapsed / count)
    total_count = g_count + y_count + r_count + b_count + o_count
    match_rate = float(total_count) / count
    print "g:{0} y:{1} r:{2} b:{3} o:{4} = {5}, out of {6}={7}".format(g_count, y_count, r_count, b_count, o_count, total_count, count, match_rate)

    print "amb: {0}, ncells: {1}".format(C.ambiguity, C.ncells)
    return results, count

# vim: et sw=2
