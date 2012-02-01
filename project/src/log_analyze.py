#!/usr/bin/env python

import sys
import numpy as np
import collections

print '#############################'
stack = []
scc = 0
_raw = {}

def canonical(name):
    return filter(lambda c: 'a' <= c <= 'z', name)
for log in sys.argv[1:]:
    print log, '->', canonical(log)
    matchdir = None
    with open(log) as f:
        for line in f:
            if 'matchdir' in line:
                matchdir = '='.join(line.strip().split('=')[1:])
            if 'Waiting for background jobs to finish' in line:
                try:
                    keys = {}
                    values = {}
                    params = eval(stack.pop())
                    keys.update(params)
                    added_error = eval(stack.pop())
                    keys.update(added_error)
                    amb, ncells = stack.pop().split(',')
                    ncells = int(ncells.split()[1])
                    keys['ncells'] = ncells
                    keys['logname'] = canonical(log)
                    keys['matchdir'] = matchdir
                    acc = stack.pop()
                    values['acc'] = float(acc.split('=')[-1])
                    total, avg = stack.pop().split(',')
                    time = float(total.split(':')[1])
#            values['t_total'] = time
                    time = float(avg.split(':')[1])
                    values['t_avg'] = time
                    stack = []
                    _k = tuple(sorted(filter(lambda (k,v): v, keys.items())))
                    if _k in _raw:
                        scc += 1
                    _raw[_k] = (keys, values)
#                    print keys,values
                except:
                    pass
            else:
                stack.append(line)
data = _raw.values()

if scc > 0:
    print '[WARN] superceded %d entries' % scc
    print '---'


n = []
def doselect(select, stream):
    collect = collections.defaultdict(list)
    print >>stream, "SELECT", select
    print >>stream, '-----------------------------'
    count = 0
    for k,v in data:
        ok = True
        for kr,vr in select.items():
            if k.get(kr) != vr:
                ok = False
        if ok:
            print >>stream, "MATCH", k,v
            for _k,_v in v.items():
                collect[_k].append(_v)
            count += 1
    print >>stream, '-----------------------------'
    print >>stream, "SUMMARY", collect
    print >>stream, '-----------------------------'
    for k, v in collect.items():
        print >>stream, k, "=", np.mean(v)
    print >>stream, "count", count
    if collect['acc']:
        n.append(len(collect['acc']))
#    print collect['acc']
#    return np.var(map(lambda x: 100*x, collect['acc']))
    return np.mean(collect['acc'])

UNIFORM = {
    'dist': 'uniform',
    'exponent': 0.5,
    'args': (0, 50**2),
}
GAUSS = {
    'dist': 'gauss',
    'exponent': 1,
    'args': (0, 23),
}
select = {
    'trees': 1,
    'amb_padding': 25,
    'vote_method': 'top_n',
}
select.update(UNIFORM)

#for ncells in [1,2,4,8,15,999]:
#    if ncells == 15:
#        print '---'
#        print 'big:',
#    elif ncells == 999:
#        print 'amb:',
#    else:
#        print '%d:' % ncells,
#    for nn in [1,2,4,8,16]:
#        select['num_neighbors'] = nn
#        select['ncells'] = ncells
#        print "%.3f" % doselect(select, open('/dev/null', 'w')),
#    print

#TEMPLATE = '/media/DATAPART2/Research/results/q5-test/matchescells(earthmine-fa10.1-culled,r=%s,d=%s),q5-test,kdtree1,threshold=70k,searchparam=1024,%stop_n'
TEMPLATE = '/media/DATAPART2/Research/results/oak-test/matchescells(oak,r=%s,d=%s),oak-test,kdtree1,threshold=70k,searchparam=1024,%stop_n'
TEMPLATE236 = '/media/DATAPART2/Research/results/oak-test/matchescells(cells),oak-test,kdtree1,threshold=70k,searchparam=1024,%stop_n'

def foo(r, nns):
    if r == 236.6:
        return TEMPLATE236 % (nns)
    else:
        return TEMPLATE % (r, r, nns)

for r in [150,236.6,350,500]:
    for nn in [1]:
        for dist, exponent, args, _ in [
            ('gauss', 1, (0,0), 25),
            ('uniform', 0.5, (0, 50**2), 75),
            ('uniform', 0.5, (0, 125**2), 150),
            ('uniform', 0.5, (0, 225**2), 250),
            ('uniform', 0.5, (0, 325**2), 350),
            ]:
            if nn > 1:
                nns = 'nn=%d,' % nn
            else:
                nns = ''
            select = {
                'args': args,
                'matchdir': foo(r, nns),
                'exponent': exponent,
                'dist': dist,
                'num_neighbors': nn,
                'ncells': 999,
            }
            x = doselect(select, open('/dev/null', 'w'))
            print r, nn, dist, args, x
        print

print "nsamples:", n
