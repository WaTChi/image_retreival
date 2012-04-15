from android import AndroidReader
import numpy as np
from info import distance

f = open('/media/DATAPART2/query5horizontal/gtLatLonNormYaw.txt')
m = {}
for line in f:
    data = line.split()
    name = data[0]
    lat, lon = float(data[1]), float(data[2])
    m[name] = (lat, lon)

dss = []
for a in AndroidReader('/media/DATAPART2/query5horizontal'):
    if a.name in m:
        d = m[a.name]
        ds = distance(a.lat, a.lon, d[0], d[1])
        ds = float(ds)
        dss.append(ds)
        dss.append(-ds)
        print ds
    else:
        print 'skipped', a.name, a.lat, a.lon
print 'mean', np.mean(dss)
print 'var', np.var(dss)
