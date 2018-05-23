#!/usr/bin/env python
# -*- coding: latin-1 -*-

from scipy import arange
from numpy import *
import time
from matplotlib.pyplot import *
from scipy.integrate import *

from modelEquations import *
from MFBfunctions import *

print "Setting up system"
### Simulation time
ti, tf = 0.0, 10e-3
tstep = 1e-5
t_eval = arange(0.0, tf, tstep)

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
#modelInput = "[0:5,0:5,0:5]"

### MFB bounding box
bb = [40, 20, 10]
boundingBox = "[0:" + str(bb[0]) + ",0:" + str(bb[1]) + ",0:" + str(bb[2]) + "]"

### Get all the compartments as
### {'i-j-k': [i,j,k,lenth,width,height]}
cmpts = compartments(modelInput)
#for k in sorted(cmpts.iterkeys()): print "%s: %s" % (k, cmpts[k])

### Check if no compartment have overlapping volumes and
### there are no gaps in the model
print "cheching geometry"
if checkGeometry(boundingBox, cmpts):
    print 'Go ahead! All good with the geometry. Take a look at it yourself'
    print 'Number of compartments:', len(cmpts)
    #plotCompartments(cmpts, bb)
else:
    exit()


### List of model objects for each compartment
cModels = []
for cname, cdim in [[k, cmpts[k]] for k in sorted(cmpts.iterkeys())]: # sorted by name
    cModels.append(mfb({'Ca':[], 'PMCA': [], 'calbindin': []},
                  name = cname,
                  dim = cdim,
                  nbrs = getNeighbours({cname: cdim}, cmpts))
                 )

c0 = '0-0-0'
cModels[0] = mfb({'Ca':[], 'HH':[], 'PMCA': [], 'calbindin': []},
                 name = c0,
                 dim = cmpts[c0],
                 nbrs = getNeighbours({c0: cmpts[c0]}, cmpts))

### print model details
#for cm in cModels: print cm.name, cm.dim, cm.nbrs, cm.nVar

### Make a list of initial index of each compartment
cmpi = initialIndex(cModels)


### Putting all compartments together
def dXdt(t, X):
    dX = []
    j=0
    for cm in cModels:
        # All compartments have V value of (0,0,0)th compartment
        cm.V = cModels[0].V

        ## Collect dX values from each compartment
        #print X
        #print X[j:j+cm.nVar]
        dX += cm.dXdt(X[j:j+cm.nVar], t)

        # Calcium Flux
        caFlux = 0
        flux = 1e4
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
print "solving the ODE"
sol = solve_ivp(dXdt, [ti, tf], X0, t_eval=t_eval)
#sol = odeint(dXdt, X0, t_eval).T
timef = time.time()
print '\nDONE! | Time elapsed: ', timef-timei

'''
### Plot some results
y = sol.y
fig, ax = subplots()
for s in y[:4]:
    figY = plot(t_eval*1e3, s, lw=1)

show()
#'''
