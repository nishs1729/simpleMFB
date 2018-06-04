#!/usr/bin/env python
from scipy import arange
from numpy import *
import time
from matplotlib.pyplot import *
from modelEquations import *
from MFBfunctions import *
from parameters import *

### Command line arguments
commandArg(sys.argv)

print "Setting up system..."
simName = 'trial/'

### Geometrical arrangement of all the compartments
### a:n   = a, a+1, a+2,...,n-1
### a:n:i = a, a+i, a+2*i,..., until (n-1)
modelInput = '''[0:2:2, 0:20:4, 0:2:3]
                [38:40:2, 0:20:4, 0:2:3]
                [0:40:5, 0:20:5, 7:10:3]
                [0:40:4, 0:20:5, 5:7:2]
                [0:40:2, 0:20:4, 3:5:2]
                [2:38:2, 0:20:2, 2:3]
                [2:38, 0:20, 0:2]
                '''

modelInput = "[0:1,0:1,0:1]"
#modelInput = "[0:2,0:1,0:1]"

### MFB bounding box
#bb = [40, 20, 10]
bb = [1]*3 #+ [1]*2
boundingBox = "[0:" + str(bb[0]) + ",0:" + str(bb[1]) + ",0:" + str(bb[2]) + "]"

### Get all the compartments as
### {'i-j-k': [i,j,k,lenth,width,height]}
cmpts = compartments(modelInput)
#for k in sorted(cmpts.iterkeys()): print "%s: %s" % (k, cmpts[k])

### Check if no compartment have overlapping volumes and
### there are no gaps in the model
print "checking geometry..."
if checkGeometry(boundingBox, cmpts):
    print 'Go ahead! All good with the geometry'
    print 'Number of compartments:', len(cmpts)
    if cmdArg['geo']:
        plotCompartments(cmpts, bb)
else:
    exit()

### List of model objects for each compartment
cModels = []
for cname, cdim in [[k, cmpts[k]] for k in sorted(cmpts.iterkeys())]: # sorted by name
    cModels.append(mfb({'Ca':[1e-7], 'PMCA': [], 'calbindin': [], 'caSensor': []},
                   name = cname,
                   dim = cdim,
                   nbrs = getNeighbours({cname: cdim}, cmpts))
                   )
'''
c0 = '0-0-0'
cModels[0] = mfb({'Ca':[], 'PMCA': [], 'calbindin': []},
                 name = c0,
                 dim = cmpts[c0],
                 nbrs = getNeighbours({c0: cmpts[c0]}, cmpts))
'''
result.solve(cModels, cmdArg=cmdArg, simName = 'trial/', flux=1e4)

### Plot some results
if cmdArg['fig']:
    fig, ax = subplots()
    for cname, c in result.data.items():
        for vname, v in result.data[cname].items():
            #if 'Ca' in vname: print vname, v[-1]
            #if vname == 'Ca':
            plot(result.t*1e3, v, lw=1, label=vname)
    legend()
    show()#'''
