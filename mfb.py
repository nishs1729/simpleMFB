#!/usr/bin/env python
# -*- coding: latin-1 -*-

from scipy import arange
from numpy import *
import time
import os
from matplotlib.pyplot import *
from scipy.integrate import *

from modelEquations import *
from MFBfunctions import *

print "Setting up system..."
### Simulation time
ti, tf = 0.0, 100e-3
tstep = 1e-4
t_eval = arange(0.0, tf, tstep)
#tcp = 100e-3

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
#modelInput = '''[0:3,0,0:3]
#                [0:3,1:3,0]
#                [0,1:3,0:3]
#                [1:2:2,1:2:2,1:2:2]'''
modelInput = "[0:5,0:5,0:5]"
modelInput = "[0:1,0:1,0:1]"

### MFB bounding box
bb = [40, 20, 10]
bb = [1]*3
boundingBox = "[0:" + str(bb[0]) + ",0:" + str(bb[1]) + ",0:" + str(bb[2]) + "]"

### Get all the compartments as
### {'i-j-k': [i,j,k,lenth,width,height]}
cmpts = compartments(modelInput)
#for k in sorted(cmpts.iterkeys()): print "%s: %s" % (k, cmpts[k])

### Check if no compartment have overlapping volumes and
### there are no gaps in the model
print "cheching geometry..."
if checkGeometry(boundingBox, cmpts):
    print 'Go ahead! All good with the geometry'
    print 'Number of compartments:', len(cmpts)
    #plotCompartments(cmpts, bb)
else:
    exit()


### List of model objects for each compartment
cModels = []
f = 1.0
for cname, cdim in [[k, cmpts[k]] for k in sorted(cmpts.iterkeys())]: # sorted by name
    cModels.append(mfb({'Ca':[f*7e-7], 'PMCA': [], 'calbindin': []},# 'PMCA': [], 'calbindin': []},
                  name = cname,
                  dim = cdim,
                  nbrs = getNeighbours({cname: cdim}, cmpts))
                 )
    f += 1.1
'''
c0 = '0-0-0'
cModels[0] = mfb({'Ca':[], 'PMCA': [], 'calbindin': []},
                 name = c0,
                 dim = cmpts[c0],
                 nbrs = getNeighbours({c0: cmpts[c0]}, cmpts))
'''

### Make a list of initial index of each compartment
cmpi = initialIndex(cModels)

'''
### print model details
for cm in cModels:
    print cm.name, cm.dim, cm.nbrs, cmpi[cm.name]
    for n in cm.nbrs:
        print cmpi[n]
#'''
tt=0.0
flux = 1e4
### Putting all compartments together
def dXdt(t, X):
    global tt
    if t>tt:
        print 't =', tt
        tt += 0.01


    dX = []
    j=0
    for cm in cModels:
        # All compartments have V value of (0,0,0)th compartment
        #cm.V = cModels[0].V

        ## Collect dX values from each compartment
        dX += cm.dXdt(X[j:j+cm.nVar], t)

        # Calcium Flux
        caFlux = 0
        for nbr in cm.nbrs:
            caFlux += flux*(X[cmpi[nbr]] - X[j]) #ADΔ(Ca)/Δx
        dX[j] += caFlux

        j += cm.nVar # increment counter to 1st element of next compartment

    return dX

### initial values
X0 = []
for cm in cModels:
    X0 += cm.X0
#print 'X0:', X0,
print 'Total number of equations:', len(X0)

### Solve ODE
timei = time.time()
print "solving the ODE..."
sol = solve_ivp(dXdt, [ti, tf], X0, t_eval=t_eval, rtol=1e-6, atol=1e-10)
#sol = odeint(dXdt, X0, t_eval).T
timef = time.time()
print '\nDONE! | Time elapsed: ', timef-timei

t = sol.t
y = sol.y

for cname, idx in cmpi.items():
    for vname, v in result.data[cname].items():
        result.data[cname][vname] = y[idx+v]

'''
### Save all the variables in files
simName = 'trial/'
for cname, c in result.data.items():
    header = 't\t\t'
    vname = [k for k in sorted(c.iterkeys())]
    for vn in vname:
        header += vn + '\t\t'
    v = array([t] + [c[k] for k in sorted(c.iterkeys())]).T
    dir = 'data/'+simName
    if not os.path.exists(dir):
        os.makedirs(dir)
    savetxt(dir+cname+'.txt', v, header=header, fmt='%.4e', delimiter='\t')
#'''


### Plot some results
if len(sys.argv) > 1:
    fig, ax = subplots()
    for cname, c in result.data.items():
        for vname, v in result.data[cname].items():
            if vname == 'Ca':
                figY = plot(t_eval*1e3, v, lw=1, label=vname)

    legend()
    show()
#'''
