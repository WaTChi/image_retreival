# The query system.
# Collection of functions that act on contexts.

import os
import geom
import math
import bisect
import sys
import time
import os.path
import query
from config import INFO
import collections
import render_tags
import posit
import pnp
import subprocess
import multiprocessing
import homographyDecomposition
import computePose
import computePose2

def get_free_mem_gb():
  txt = subprocess.Popen(['free', '-g'], stdout=subprocess.PIPE).communicate()[0]
  return int(txt.split('\n')[2].split()[-1])

# based on memory usage
# rule of thumb:
#   we need 1.3gb of free memory per thread
#   and should allow 10gb for disk cache
t_avg = multiprocessing.cpu_count()
tmax = multiprocessing.cpu_count()
last_sampled = 0
thread_mem_gb = 2.0
def estimate_threads_avail():
  if os.getenv('NUM_THREADS'):
    return int(os.getenv('NUM_THREADS'))
  global t_avg, last_sampled
  if time.time() - last_sampled > 5:
    last_sampled = time.time()
    gb_free = get_free_mem_gb()
    t = int((gb_free - 10.0)/thread_mem_gb)
    t_avg = t*.3 + t_avg*.7
    INFO("I think we have enough memory for %d threads" % t_avg)
  t = max(1, min(tmax, int(t_avg)))
  if os.getenv('MAX_THREADS'):
    return min(t, int(os.getenv('MAX_THREADS')))
  return t

from config import *
import info
import numpy as np
import corr
import query1GroundTruth
import query2Groundtruth
import query5horizGroundTruth
import query5vertGroundTruth
import oakland1GroundTruth
import groundtruthB
import groundtruthG
import groundtruthO
import groundtruthR
import groundtruthY
import query4GroundTruth
import cory25GroundTruth
import emeryvilleGroundTruth
import util
import util_cs188

from multiprocessing import Pool, cpu_count

class MultiprocessExecution:
  pool = None

  def __enter__(self):
    MultiprocessExecution.pool = Pool(multiprocessing.cpu_count())

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

def derive_key(C, closest_cells, name):
    return (name,C.params['num_neighbors'],C.params['trees']) + tuple(sorted(map(lambda (cell, dist): cell, closest_cells)))

class DiscreteCell(object):
  def __init__(self):
    self.count = 0
    self.images = set()

  def incr(self):
    self.count += 1

  def register(self, i):
    self.images.add(i)

  def __repr__(self):
    return "%d over %d" % (self.count, len(self.images))

def weight_by_coverage(C, Q, matches):
  area = collections.defaultdict(DiscreteCell)
  scores = collections.defaultdict(int)

  def discretize(lat, lon):
    EARTH_DIAMETER = 12756000 # meters
    SIZE = 3.0 # meters


    size = 360.0*SIZE/EARTH_DIAMETER
    return size*int(lat/size), size*int(lon/size)

  pm = C.pixelmap

  # put all 3d points of db image into map
  # for possible matching
  preregister = False
  if preregister:
    for image, _, _ in matches:
      locs = pm.open(image + 'sift.txt')
      for loc in locs.values():
        if loc:
          lld = discretize(loc['lat'], loc['lon'])
          area[lld].register(image)

  for target in matches:
    image, score, pairs = target

    for p in pairs:
      x, y = p['db'][:2]
      loc = pm.get(image + 'sift.txt', x, y)
      if loc is None:
        # increment individual score by 1
        scores[image] += 1
      else:
        # increment area score for loc by 1
        lld = discretize(loc['lat'], loc['lon'])
        area[lld].incr()
        area[lld].register(image)

  print "Area map len", len(area)
  # tally area contributions
  for k, cell in area.iteritems():
    for i in cell.images:
      scores[i] += cell.count

  return distance_sort(C, Q, scores.items())

def weight_by_distance(C, Q, matches):

  def weighting_function(dist):
    # exponential decay with w=0.36 @ 20 meters
    ONE_OVER_E_POINT = 20.0
    return math.e**(-dist/ONE_OVER_E_POINT)

  def contrib_function(contributer):
    # squared contribution
    return contributer[1]**2

  def contribution(target, contributer):
    lat, lon = getlatlonfromdbimagename(C, target[0])
    lat2, lon2 = getlatlonfromdbimagename(C, contributer[0])
    dist = info.distance(lat, lon, lat2, lon2)
    return contrib_function(contributer) * weighting_function(dist)
  newmatches = []
  for target in matches:
    score = 0.0
    for contributer in matches:
      score += contribution(target, contributer)
    newmatches.append((target[0], score))
  print "BEFORE", distance_sort(C, Q, matches)[0]
  print "AFTER", distance_sort(C, Q, newmatches)[0]
  return distance_sort(C, Q, newmatches)

def amb_cutoff(C, Q, matches):
  def admit(match):
    line = match[0][:-8]
    lat, lon = getlatlonfromdbimagename(C, line)
    dist = info.distance(lat, lon, Q.query_lat, Q.query_lon)
    if dist < C.amb_cutoff + C.amb_padding:
      return True
    return False
  return filter(admit, matches)

def distance_sort(C, Q, matches):
  def extract_key(match):
    line = match[0]
    lat, lon = getlatlonfromdbimagename(C, line)
    return (match[1], -info.distance(lat, lon, Q.query_lat, Q.query_lon))
  return sorted(matches, key=extract_key, reverse=True)

cache = {}
def match(C, Q):
    if C.shuffle_cells:
      C._dbdir = None
    if C.override_cells:
      INFO('override cells')
      cells_in_range = [(c,0) for c in C.override_cells]
    else:
      # compute closest cells
      closest_cells = util.getclosestcells(Q.query_lat, Q.query_lon, C.dbdir)
      if C.restrict_cells:
        closest_cells = filter(lambda c: c[0] in C.restrict_cells, closest_cells)
      cells_in_range = [(cell, dist)
        for cell, dist in closest_cells[0:C.ncells]
          if dist < C.cellradius + C.ambiguity + C.matchdistance]
    INFO('Using %d cells' % len(cells_in_range))
    if C.shuffle_cells:
      import reader
      sr = reader.get_reader('sift')
      supercell = sr.get_supercelldir(
        C.dbdir,
        [c for (c,d) in cells_in_range],
        C.overlap_method)
      C._dbdir = supercell

    if not cells_in_range:
        raise LocationOutOfRangeError

    # cache for fuzz runs
    if C.cacheEnable:
        key = derive_key(C, cells_in_range, Q.siftname)
        if key in cache:
            print 'cache hit'
            return cache[key]
        else:
            print 'cache miss'

    # compute output file paths for the cells

    cellpath = [c for c,d in cells_in_range]
    if C.one_big_cell:
      INFO('Using 1 big cell (%d union)' % len(cells_in_range))
      outputFilePaths = [os.path.join(C.matchdir, Q.siftname + ',' + getcellid(cellpath) + ".res")]
      cellpath = [cellpath]
    else:
      outputFilePaths = []
      for cell, dist in cells_in_range:
          if ',' in cell:
            latcell, loncell = cell.split(',')
            latcell = float(latcell)
            loncell = float(loncell)
          else:
            latcell, loncell = 0,0
          actualdist = info.distance(Q.query_lat, Q.query_lon, latcell, loncell)
          outputFilePath = os.path.join(C.matchdir, Q.siftname + ',' + cell + ',' + str(actualdist)  + ".res")
          outputFilePaths.append(outputFilePath)

    # start query
    query.run_parallel(C, Q, cellpath, outputFilePaths, estimate_threads_avail())

    # combine results
    if C.spatial_comb:
      comb_matches = corr.combine_spatial(outputFilePaths)
    else:
      comb_matches = corr.combine_matches(outputFilePaths)

    #geometric consistency reranking
    if C.disable_filter_step:
      imm = condense2(sorted(comb_matches.iteritems(), key=lambda x: len(x[1]), reverse=True))
      rsc_ok = True
    else:
      imm, rsc_ok = rerank_ransac(comb_matches, C, Q)

    if C.weight_by_coverage:
      ranked = weight_by_coverage(C, Q, imm)
    elif C.weight_by_distance:
      ranked = weight_by_distance(C, Q, imm)
    else:
      ranked = distance_sort(C, Q, imm)

    # top 1
    stats = check_topn_img(C, Q, ranked, 1)

    # return statistics and top result
    matchedimg = ranked[0][0]
    matches = comb_matches[matchedimg + 'sift.txt']
    if C.cacheEnable:
        cache[key] = (stats, matchedimg, matches, ranked)
    if C.match_callback:
        C.match_callback(C, Q, stats, matchedimg, ranked, cells_in_range, rsc_ok)

    # compute homography and draw images maybe
    if MultiprocessExecution.pool:
      MultiprocessExecution.pool.apply_async(compute_hom, [C.pickleable(), Q, ranked, comb_matches])
    else:
      compute_hom(C, Q, ranked, comb_matches)

    ### Query Pose Estimation ###
    match = any(check_img(C, Q, ranked[0]))
    if (C.solve_pose and match and Q.name not in C.pose_remove) or C.pose_param['solve_bad']:
        #computePose.draw_dbimage(C, Q, matchedimg, match)
        if MultiprocessExecution.pool:
            MultiprocessExecution.pool.apply_async(computePose.estimate_pose, [C.pickleable(), Q, matchedimg, match])
        else:
            computePose.estimate_pose(C, Q, matchedimg, match)

    # done
    return stats, matchedimg, matches, ranked

def getlatlonfromdbimagename(C, dbimg):
    if C.QUERY == 'emeryville' or C.QUERY == 'cory-25' or C.QUERY == 'cory-2' or C.QUERY == 'cory-5':
        return 0,0
    clat = float(dbimg.split(",")[0])
    clon = float(dbimg.split(",")[1][0:-5])
    return clat, clon
    
# top entry is ranked_matches[0], etc
def compute_hom(C, Q, ranked_matches, comb_matches):
    match1 = any(check_img(C, Q, ranked_matches[0]))
    if not C.compute_hom:
      if match1:
        return
      if not C.log_failures:
        return
    if C.put_into_dirs:
        udir = os.path.join(C.resultsdir, Q.name)
    else:
        udir = C.resultsdir
    if not os.path.exists(udir):
        os.makedirs(udir)
    i = 0
    data = {}
    for matchedimg, score, pairs in ranked_matches[:C.max_matches_to_analyze]:
        i += 1
        if C.stop_on_homTrue and data.get('success'):
            break # we are done (found good homography)
        clat, clon = getlatlonfromdbimagename(C, matchedimg)
        matchimgpath = None
        # XXX this sucks, since we don't have a db image abstraction
        for ext in ['.jpg', '.JPG', '.png', '.PNG']:
          p = os.path.join(C.dbdump, '%s%s' % (matchedimg, ext))
          if os.path.exists(p):
            matchimgpath = p
        assert matchimgpath
        match = any(check_img(C, Q, ranked_matches[i-1]))

#        matches = db_matches
        # rematch for precise fit
        db_matches = comb_matches[matchedimg + 'sift.txt']
        matchsiftpath = os.path.join(C.dbdump, matchedimg + 'sift.txt')
        matches = corr.rematch(C, Q, matchsiftpath)
#        matches1 = matches
#        rp_matches = corr.rematch(C, Q, matchsiftpath)

        # concat db matches
        matches.extend(db_matches)

        # find homography
        rsc_matches, H, inliers = corr.find_corr(matches, hom=True, ransac_pass=True, data=data)
        rsc_inliers = np.compress(inliers, rsc_matches).tolist()
        u = corr.count_unique_matches(rsc_inliers)

        if C.drawtopcorr or (not match and C.log_failures):
          # draw picture
          matchoutpath = os.path.join(udir, Q.name + ';match' + str(i) + ';gt' + str(match)  + ';hom' + str(data.get('success')) + ';uniq=' + str(u) + ';inliers=' + str(float(sum(inliers))/len(matches)) + ';' + matchedimg + '.jpg')
#          try:
          corr.draw_matches(C, Q, matches, rsc_matches, H, inliers, matchimgpath, matchoutpath, matchsiftpath, C.show_feature_pairs)
#          except IOError, e:
#            INFO(e)

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

          
def condense2(list):
    return map(lambda x: (x[0][:-8], len(x[1]), x[1]), list)

def condense3(list):
    return map(lambda x: (x[1][:-8], x[0], x[2]), list)

def rerank_ransac(counts, C, Q):
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)

    if C.amb_cutoff:
      assert C.amb_cutoff > 1
      old = sorted_counts
      sorted_counts = amb_cutoff(C, Q, sorted_counts)
      print "amb cutoff: %d -> %d" % (len(old), len(sorted_counts))
      if not sorted_counts:
        sorted_counts = old

    # tuples of (len, siftfile, matches) ordered [(1..), (2..), (3..)]
    reranked = []
    num_filt = 0

    for siftfile, matches in sorted_counts:
      # the bottom g items in 'reranked' are in their final order
      g = len(reranked) - bisect.bisect(reranked, (len(matches), None, None))
      if g >= C.ranking_min_consistent or num_filt >= C.ranking_max_considered:
        if C.verbosity > 0:
          INFO('stopped after filtering %d' % num_filt)
        break
      num_filt += 1

      # perform ransac
      F, inliers = corr.find_corr(matches)
      if any(inliers):
        m = np.compress(inliers, matches)
        bisect.insort(reranked, (len(m), siftfile, m))

    if not reranked:
      INFO('W: no db matches passed ransac, not filtering')
      return condense2(sorted_counts), False
    else:
      reranked.reverse()
      return condense3(reranked), True

def combine_vote(counts):
    sorted_counts = sorted(counts.iteritems(), key=lambda x: len(x[1]), reverse=True)
    return condense2(sorted_counts), False

def check_img(C, Q, entry):
    g,y,r,b,o = 0,0,0,0,0
    if C.QUERY == 'query1' or C.QUERY == 'query1-m' or C.QUERY == 'query1-a':
        g += check_truth(Q.name, entry[0], query1GroundTruth.matches)
    elif C.QUERY == 'oakland1' or C.QUERY == 'oak-test':
        g += check_truth(Q.name, entry[0], oakland1GroundTruth.matches)
    elif C.QUERY == 'query3' or C.QUERY == 'query3a':
        g += check_truth(Q.name, entry[0], groundtruthG.matches)
        y += check_truth(Q.name, entry[0], groundtruthY.matches)
        r += check_truth(Q.name, entry[0], groundtruthR.matches)
        b += check_truth(Q.name, entry[0], groundtruthB.matches)
        o += check_truth(Q.name, entry[0], groundtruthO.matches)
    elif C.QUERY == 'query2':
        g += check_truth(Q.name, entry[0], query2Groundtruth.matches)
    elif C.QUERY == 'query4' or C.QUERY == 'query4-cropped' or C.QUERY == 'query4a' or C.QUERY == 'q4-test':
        g += check_truth(Q.name, entry[0], query4GroundTruth.matches)
    elif C.QUERY == 'query5horizontal' or C.QUERY == 'q5-test':
        g += check_truth(Q.name, entry[0], query5horizGroundTruth.matches)
    elif C.QUERY == 'query5vertical':
        g += check_truth(Q.name, entry[0], query5vertGroundTruth.matches)
    elif C.QUERY == 'cory-25' or C.QUERY == 'cory-2' or C.QUERY == 'cory-5':
        g += check_truth(Q.name, entry[0], cory25GroundTruth.matches)
    elif C.QUERY == 'emeryville':
        g += check_truth(Q.name, entry[0], emeryvilleGroundTruth.matches)
    else:
        return [0,0,0,0,0]
    return [g > 0, y > 0, r > 0, b > 0, o > 0]

def check_topn_img(C, Q, dupCountLst, topnres=1):
    record = [0]*5
    for entry in dupCountLst[0:topnres]:
        new = check_img(C, Q, entry)
        record = map(lambda a,b: a + b, record, new)
    return map(bool, record)

def dump_combined_matches(C, Q, stats, matchedimg, matches, cells_in_range, rsc_ok):
    # For Aaron's analysis
    table = {}
    for line in open(os.path.join(C.dbdir, 'cellmap.txt')):
        a, b = line.split()
        table[b] = int(a)
    def cellsetstr(cells):
        cells = sorted(map(lambda (cell, dist): str(table[cell]), cells))
        return '-'.join(cells)
    outputFilePath = os.path.join(C.matchdir, C.aarondir, Q.siftname + ',combined,' + cellsetstr(cells_in_range) + ".res")
    d = os.path.dirname(outputFilePath)
    if not os.path.exists(d):
        os.makedirs(d)
    def save(outputFilePath):
        with open(outputFilePath, 'w') as outfile:
            if rsc_ok:
              outfile.write('ransac_ok\n')
            else:
              outfile.write('ransac_failed\n')
            for matchedimg, score in matches:
                outfile.write(str(score))
                outfile.write('\t')
                outfile.write(matchedimg)
                outfile.write('\n')
    save_atomic(save, outputFilePath)

def characterize(C):
    print >>sys.stderr,"start"
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
    try:
      open(C.pose_param['pose_file'],'w').close()
      open(C.pose_param['extras_file'],'w').close()
    except Exception, e:
      INFO(e)
    for Q in C.iter_queries():
        try:
          Q.datafile = os.path.join(C.pose_param['resultsdir'],'data_'+Q.name+'.txt')
          open(Q.datafile,'w').close()
        except Exception, e:
          INFO(e)
        if C.verbosity>0:
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
    print C.added_error
    print C.params
    return results, count

# vim: et sw=2
