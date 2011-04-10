#!/usr/bin/env python

import commands
from sys import stdin

bin = '../bins/' + commands.getoutput('uname -m') + '-' + commands.getoutput('uname -s')

spots = []
count = []
data  = []

def angle(x,y):
        """angular separation between two pixels on a HealPix map (in arcmin)"""
        
        exe = bin + '/pxltool -a ../maps/white.fits'
        out = commands.getoutput('%s %i %i' % (exe,x,y))
        
        return float(out)


# find cold spot clusters
for l in stdin:
        if l[0] == '#': continue
        token = l.split(); w = float(token[0])/2.0; p = int(token[1])
        
        k = 0
        
        for s in spots:
                if p in s: break
                if min([angle(i,p) for i in s]) < w:
                        s.append(p)
                        break
                k += 1
        else:
                spots.append([p])
                count.append(0)
        count[k] += 1
        
        data.append((k,l))

# rank the coldspots by filters triggered
idx = range(len(spots)); idx.sort(lambda a,b: cmp(count[b],count[a]))
rank = range(len(spots)); rank.sort(lambda a,b: cmp(idx[a],idx[b]))

# output data
for k,l in data:
        print l.strip(), rank[k]+1

print '# cold spot clusters:', [spots[i] for i in idx]
