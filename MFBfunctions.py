from matplotlib.pyplot import *
from scipy.integrate import *
from numpy import *
import os
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from parameters import *

### Takes in location and size of each compartment
### and returns a dictionary of compartment name and
### corresponding location ad size array
def compartments(modelInput):
    modelInput = [a.strip() for a in modelInput.strip().split("\n")]
    cmpts = {}
    for a in modelInput:
        a = a.strip('[').strip(']').split(',')
        #print 'a:', a

        arr, temp = [], []
        for i in a[:3]:
            if ':' in i:
                t = [int(a) for a in i.split(':')]
                try:
                    arr.append(range(t[0],t[1],t[2]))
                    temp.append(t[2])
                except IndexError:
                    arr.append(range(t[0],t[1]))
                    temp.append(1)
            else:
                arr.append([int(i)])
                temp.append(1)

        arr += temp

        for i in arr[0]:
            for j in arr[1]:
                for k in arr[2]:
                    mName = str(i) + '-' + str(j) + '-' + str(k)
                    mVal = [i,j,k] + arr[3:]
                    cmpts.update({mName: mVal})
    return cmpts

### Check if no compartment have overlapping volumes and
### there are no gaps in the model.
### comp is the bounding box of all compartments broken into
### smallest cubic unit. For eg. "[0:3,0:3,0:3]".
### cmpts is the dictionary of all the compartments.
def checkGeometry(boundingBox, cmpts):
    sCubes = compartments(boundingBox)

    gap = []
    overlap = []
    for sc in sCubes.values():
        temp = []
        i = 0
        for c in cmpts.values():
            if (sc[0]  >=c[0]      and sc[1]  >=c[1]      and sc[2]  >=c[2] and \
                sc[0]+1<=c[0]+c[3] and sc[1]+1<=c[1]+c[4] and sc[2]+1<=c[2]+c[5]):
                i += 1
                temp.append(c)

        if i>=2:
            overlap.append([sc[:3], temp])
        if i==0:
            gap.append(sc[:3])

    if gap==[] and overlap==[]:
        return True
    else:
        print 'Error in compartment geometry!'
        print 'gaps:', gap, len(gap)
        print 'overlaps:', overlap
        return False


### Check if two compartments a amd b are neighbours
### Returns False if they are not neighbours
### Returns shared surface area if they are neighbours
def isNeighbour(a,b):
    nbr = (a[0] <= b[0]+b[3] and a[0]+a[3] >= b[0]) and \
          (a[1] <= b[1]+b[4] and a[1]+a[4] >= b[1]) and \
          (a[2] <= b[2]+b[5] and a[2]+a[5] >= b[2])


    if nbr:
        if a[0] == b[0]+b[3] or a[0]+a[3] == b[0]:
            dy = min(abs(a[1]+a[4] - b[1]), abs(b[1]+b[4] - a[1]))
            dz = min(abs(a[2]+a[5] - b[2]), abs(b[2]+b[5] - a[2]))

            if dy*dz==0:
                return False
            else:
                return dy*dz

        elif a[1] == b[1]+b[4] or a[1]+a[4] == b[1]:
            dx = min(abs(a[0]+a[3] - b[0]), abs(b[0]+b[3] - a[0]))
            dz = min(abs(a[2]+a[5] - b[2]), abs(b[2]+b[5] - a[2]))

            if dx*dz==0:
                return False
            else:
                return dx*dz

        elif a[2] == b[2]+b[5] or a[2]+a[5] == b[2]:
            dx = min(abs(a[0]+a[3] - b[0]), abs(b[0]+b[3] - a[0]))
            dy = min(abs(a[1]+a[4] - b[1]), abs(b[1]+b[4] - a[1]))

            if dx*dy==0:
                return False
            else:
                return dx*dy

    else:
        return False

### Get all neighbours of a compartment c
### cmpts is the dictionary of all the compartments
def getNeighbours(c, cmpts):
    k, v = c.items()[0]
    nbrs = []
    for nk,nv in cmpts.items():
        if not k==nk:
            area = isNeighbour(v, nv)
            if area:
                nbrs.append(nk)
    return nbrs


### Get all vertices of a compartment
def getVertices(c):
    v = [c[:3]]
    v.append([c[0]+c[3], c[1],      c[2]])
    v.append([c[0]+c[3], c[1]+c[4], c[2]])
    v.append([c[0],      c[1]+c[4], c[2]])
    v.append([c[0],      c[1],      c[2]+c[5]])
    v.append([c[0]+c[3], c[1],      c[2]+c[5]])
    v.append([c[0]+c[3], c[1]+c[4], c[2]+c[5]])
    v.append([c[0],      c[1]+c[4], c[2]+c[5]])

    # list of sides' polygons of figure
    verts = [[v[0], v[1], v[2], v[3]],
             [v[4], v[5], v[6], v[7]],
             [v[0], v[1], v[5], v[4]],
             [v[2], v[3], v[7], v[6]],
             [v[1], v[2], v[6], v[5]],
             [v[4], v[7], v[3], v[0]],
             [v[2], v[3], v[7], v[6]]]

    return array(verts)

### Draw each compartment of the bouton for visual verification
def plotCompartments(cmpts, bb):
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')

    # set x,y,z limits
    ax.set_xlim(0, bb[0])
    ax.set_ylim(0, bb[1])
    ax.set_zlim(0, bb[2])

    ax.axis('off')

    # plot sides
    for c in cmpts.values():
        verts = getVertices(c)
        collection = Poly3DCollection(verts, linewidths=0.5,
                                    edgecolors='g', alpha=.4)
        collection.set_facecolor('cyan')
        ax.add_collection3d(collection)
        #axis('equal')
    show()


### Make a list of initial index of each compartment
def initialIndex(cModels):
    cmpi = {}
    i = 0
    for cm in cModels:
        cmpi.update({cm.name: i})
        i += cm.nVar
    return cmpi

### Command line arguments
def commandArg(argv):
    for a in argv[1:]:
        a = a.split('=')
        cmdArg.update({a[0]: float(a[1])})


### The solution class
class solution:
    data = {}
    t = 0

    ### Putting all compartments together
    def dXdt(self, t, X):
        if t>self.t:
            print 't =', self.t
            self.t += self.tcp/5

        dX = []
        j=0
        for cm in self.cModels:
            ## All compartments have V value of (0,0,0)th compartment
            #cm.V = cModels[0].V

            ## Collect dX values from each compartment
            dX += cm.dXdt(X[j:j+cm.nVar], t)

            ## Calcium Flux
            caFlux = 0
            for nbr in cm.nbrs:
                caFlux += self.flux*(X[cmpi[nbr]] - X[j])
            dX[j] += caFlux

            j += cm.nVar # increment counter to 1st element of next compartment

        return dX

    ### Solving the ODEs
    def solve(self, cModels, cmdArg={}, simName = 'trial/', flux=1e4):
        self.cModels = cModels
        self.flux = flux

        ## Make a list of initial index of each compartment
        cmpi = initialIndex(cModels)

        ## Simulation time
        ti, tf = 0.0, cmdArg['tf']
        tstep = cmdArg['tstep']
        self.tcp = cmdArg['tcp']

        ## initial values
        X0 = []
        for cm in cModels:
            X0 += cm.X0
        print 'Total number of equations:', len(X0)

        ## Solve ODE
        timei = time.time()
        print "solving the ODE..."

        tinterval = []
        aa = arange(ti, tf, self.tcp)
        for i in range(len(aa[:-1])):
            tinterval += [[aa[i], aa[i+1]]]
        tinterval += [[aa[-1], tf]]

        ## Saving data till current checkpoint in files
        if cmdArg['save']:
            dir = 'data/'+cmdArg['simName']
            if not os.path.exists(dir):
                os.makedirs(dir)
            fname = dir+cm.name+'.txt'
            file = open(fname, 'w')

        temp = 1
        for ti, tf in tinterval:
            t_eval = linspace(ti, tf, round((tf - ti)/tstep) + 1)[:-1]

            sol = solve_ivp(self.dXdt, [ti, tf], X0, t_eval=t_eval,
                            rtol=cmdArg['rtol'], atol=cmdArg['atol'])

            ## Last values of sol as X0
            X0 = sol.y.T[-1]

            ## Saving data till current checkpoint in files
            if cmdArg['save']:
                tsavei = time.time()
                for cm in cModels:
                    header = 't\t\t'
                    vname = sorted(self.data[cm.name].iteritems(), key=lambda (k,v): (v,k))
                    for vn in vname:
                        header += vn[0] + '\t\t'

                    v = concatenate(([sol.t], sol.y[cmpi[cm.name]:cmpi[cm.name]+cm.nVar])).T
                    if temp:
                        savetxt(file, v, header=header, fmt='%.4e', delimiter='\t')
                    else:
                        savetxt(file, v, fmt='%.4e', delimiter='\t')

                tsavef = time.time()
                print 'Time taken for saving:', tsavef - tsavei

            ## Adding sol to solution till previous checkpoint
            if temp:
                t = sol.t
                y = sol.y
                temp = 0
            else:
                t = concatenate((t, sol.t))
                y = concatenate((y, sol.y), axis=1)

        file.close()
        ## Organise data in result.data dictionary
        for cname, idx in cmpi.items():
            for vname, v in self.data[cname].items():
                self.data[cname][vname] = y[idx+v]
        self.t = t
