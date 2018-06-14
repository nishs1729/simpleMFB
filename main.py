#!/usr/bin/env python
from mfb import *

timei = time.time() # for calculating total run time

### Command line arguments
commandArg(sys.argv)

print "Setting up system..."
### Enter model description
modelDesc = modelDesc200

### MFB bounding box
bBox = boundingBox(modelDesc)

### Get all the compartments as
### {'i-j-k': [i,j,k,lenth,width,height]}
cmpts = compartments(modelDesc)
#for k, v in cmpts.items(): print k, v

### Compartment list for specific model type
cm = {
    'HH': ['0-0-0'],
    'PMCA': [], #cSurf(cmpts)
    'VDCC': [], #hexPoints(9, 500, bBox),
    'caSensor': [],
    'calbindin': [] #cmpts.keys()
}

### Check if no compartment have overlapping volumes and
### there are no gaps in the model
if cmdArg['geo']:
    print "Compartment geometry:"
    if checkGeometry(bBox, cmpts):
        print Fore.GREEN + 'Go ahead! All good with the geometry\n' + Fore.RESET
        print 'Compartments:\t', Fore.GREEN, len(cmpts), Style.RESET_ALL
        plotCompartments(cmpts, bBox, cm)
    else:
        exit()
else:
    print 'Compartments:\t', Fore.GREEN, len(cmpts), Style.RESET_ALL

### List of model objects for each compartment
cModels = getModels(cmpts, cm)

### Create a solution object
result = solution(cModels, cmpts)

### Solve the equations
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    result.solve()

### Display total time taken
timef = time.time()
print '\n\tTotal time:', Fore.RED, timef-timei, Style.RESET_ALL

### Plot some results
if cmdArg['fig']:
    fig, ax = plt.subplots()
    for cname, c in result.data.items():
        for vname, v in result.data[cname].items()[:]:
            #if 'Ca' in vname: print vname, v[-1]
            #if vname == 'Ca':
            plt.plot(result.t*1e3, v, lw=1, label=vname)

    if cmdArg['vfile']:
        v = getV()
        plt.plot(result.t*1e3, v(result.t), lw=1, label='V')

    plt.legend()
    plt.show()
#"""
