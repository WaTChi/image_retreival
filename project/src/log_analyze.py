#!/usr/bin/env python

import sys
import numpy as np
import collections

stack = []
_raw = {}
with open(sys.argv[1]) as f:
    for line in f:
        if 'Waiting for background jobs to finish' in line:
            keys = {}
            values = {}
            params = eval(stack.pop())
            keys.update(params)
            added_error = eval(stack.pop())
            keys.update(added_error)
            amb, ncells = stack.pop().split(',')
            ncells = int(ncells.split()[1])
            keys['ncells'] = ncells
            acc = stack.pop()
            values['acc'] = float(acc.split('=')[-1])
            total, avg = stack.pop().split(',')
            time = float(total.split(':')[1])
            values['t_total'] = time
            time = float(avg.split(':')[1])
            values['t_avg'] = time
            stack = []
            _raw[tuple(sorted(keys.items()))] = (keys, values)
        else:
            stack.append(line)
data = _raw.values()

select = {
    'dist': 'gauss',
    'trees': 4,
    'ncells': 4,
    'num_neighbors': 4,
}

collect = collections.defaultdict(list)

print '-----------------------------'
print "SELECT", select
print '-----------------------------'
count = 0
for k,v in data:
    ok = True
    for kr,vr in select.items():
        if k[kr] != vr:
            ok = False
    if ok:
        print "MATCH", k,v
        for _k,_v in v.items():
            collect[_k].append(_v)
        count += 1
print '-----------------------------'
print "SUMMARY", collect
print '-----------------------------'
for k, v in collect.items():
    print k, "=", np.mean(v)
