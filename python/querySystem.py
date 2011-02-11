import queryContext as context
import os.path
import time
from config import save_atomic, INFO

context.QUERY = 'query3'
context.params.update({
  'checks': 1024,
  'trees': 1,
  'vote_method': 'filter',
  'num_neighbors': 1,
  'dist_threshold': 70000,
  'confstring': '',
})

count = 0
start = time.time()
fuzzydir = os.path.join(context.maindir, 'fuzzylocs/%s' % context.QUERY)
results = {}
topnresults = [1,2,5,10]
for n in topnresults:
    results[n]=0

def match_callback(siftfile, matchdir, stats, matchedimg, matches, cells_in_range):
    # For Aaron's analysis
    table = {}
    for line in open(os.path.join(context.dbdir, 'cellmap.txt')):
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

    # compile statistics
    for n in topnresults:
        result = context.check_topn_img(siftfile, matches, n)
        results[n] = reduce(lambda x,y: x or y, result)
    global count
    count += 1
    INFO('speed is %f' % ((time.time()-start)/count))
    for n in topnresults:
        print "matched {0}\t out of {1}\t in the top {2}\t amb: {3}, ncells:{4}".format(results[n], count, n, context.ambiguity, context.ncells)

context.ambiguity = 75
context.locator_function = context.skew_location
context.match_callback = match_callback
context.cacheEnable = 1
context.characterize()

for n in topnresults:
    total_count = results[n]
    match_rate = float(total_count) / count
    print "matched {0}\t out of {1}\t = {2}\t in the top {3}".format(total_count, count, match_rate, n)
