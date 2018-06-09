#!/usr/bin/env python
from mfb import *

timei = time.time() # for calculating total run time

### Command line arguments
commandArg(sys.argv)

print "Setting up system..."
simName = 'trial/'

### Geometrical arrangement of all the compartments
### a:n   = a, a+1, a+2,...,n-1
### a:n:i = a, a+i, a+2*i,..., until (n-1)
modelDesc = '''[0:2:2, 0:20:4, 0:2:3]
                [38:40:2, 0:20:4, 0:2:3]
                [0:40:5, 0:20:5, 7:10:3]
                [0:40:4, 0:20:5, 5:7:2]
                [0:40:2, 0:20:4, 3:5:2]
                [2:38:2, 0:20:2, 2:3]
                [2:38, 0:20, 0:2]
                '''

modelDesc = "[0:1,0:1,0:1]"
#modelDesc = "[0:2,0:2,0:2]"
#modelDesc = "[0:5,0:5,0:5]"

### MFB bounding box
bb = [40, 20, 10]
bb = [1]*3 #+ [1]*2

bb = [str(a) for a in bb]
boundingBox = "[0:" + bb[0] + ",0:" + bb[1] + ",0:" + bb[2] + "]"

### Get all the compartments as
### {'i-j-k': [i,j,k,lenth,width,height]}
cmpts = compartments(modelDesc)
#for k, v in cmpts.items(): print k, v

### Compartment list for specific model type
cm = {
    'cHH': ['0-0-0'],
    'cVDCC': ['1-0-0', '1-1-0', '1-0-1', '1-1-1', '0-0-1', '0-1-0', '0-1-1'],
    'cPMCA': [],
    'cPMCASensor': []
}

### Check if no compartment have overlapping volumes and
### there are no gaps in the model
print "Compartment geometry:"
if checkGeometry(boundingBox, cmpts):
    print Fore.GREEN + 'Go ahead! All good with the geometry\n' + Fore.RESET
    print 'Compartments:\t', Fore.GREEN, len(cmpts), Style.RESET_ALL
    if cmdArg['geo']:
        plotCompartments(cmpts, bb, cm)
else:
    exit()


### List of model objects for each compartment
cModels = getModels(cmpts, cm)
#print cModels

### Create a solution object
result = solution(cModels, cmpts, cmdArg, simName='trial/')

### Solve the equations
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    result.solve()



timef = time.time()
print '\n\tTotal time:', Fore.RED, timef-timei, Style.RESET_ALL

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
#"""
