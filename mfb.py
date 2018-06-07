#!/usr/bin/env python
from modelEquations import *
from MFBfunctions import *
from parameters import *
from solution import *

timei = time.time()
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

#modelInput = "[0:1,0:1,0:1]"
modelInput = "[0:2,0:2,0:2]"

### MFB bounding box
#bb = [40, 20, 10]
bb = [2]*3 #+ [1]*2
boundingBox = "[0:" + str(bb[0]) + ",0:" + str(bb[1]) + ",0:" + str(bb[2]) + "]"

### Get all the compartments as
### {'i-j-k': [i,j,k,lenth,width,height]}
cmpts = compartments(modelInput)
#for k, v in cmpts.items(): print k, v

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
cModels = od()
for cname, cdim in cmpts.items():
    cModels.update({cname: mfb({'Ca':[1e-7], 'PMCA': [], 'HH': []},
                               name = cname,
                               dim = cdim)
                  })

'''
c0 = '0-0-0'
cModels.update({c0: mfb({'Ca':[100e-7], 'PMCA': [], 'calbindin': []},
                name = c0,
                dim = cmpts[c0],
                nbrs = getNeighbours(c0, cmpts[c0], cmpts))}
              )
'''

### Create a solution object
result = solution(cModels, cmpts, cmdArg, simName='trial/')

### Solve the equations
result.solve()



timef = time.time()
print '\n\tTotal time:', timef-timei

### Plot some results
if cmdArg['fig']:
    fig, ax = plt.subplots()
    for cname, c in result.data.items():
        for vname, v in result.data[cname].items()[:5]:
            #if 'Ca' in vname: print vname, v[-1]
            #if vname == 'Ca':
            plt.plot(result.t*1e3, v, lw=1, label=vname)
    plt.legend()
    plt.show()
