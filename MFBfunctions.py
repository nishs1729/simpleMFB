from parameters import *
from modelEquations import *

### Takes in location and size of each compartment
### and returns a dictionary of compartment name and
### corresponding location ad size array
def compartments(modelInput):
    modelInput = [a.strip() for a in modelInput.strip().split("\n")]
    cmpts = od()
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
#@jit
def isNeighbour(a,b):
    nbr = (a[0] <= b[0]+b[3] and a[0]+a[3] >= b[0]) and \
          (a[1] <= b[1]+b[4] and a[1]+a[4] >= b[1]) and \
          (a[2] <= b[2]+b[5] and a[2]+a[5] >= b[2])

    if nbr:
        ## distance between the compartments, d
        mida = a[:3] + a[3:]/2.0
        midb = b[:3] + b[3:]/2.0
        #d = np.sqrt(np.sum((mida-midb)**2))
        d = np.sqrt(reduce(op.add, map(lambda x, y: (x-y)**2, mida, midb)))

        ## If a and b have an overlapping surface, return area and d
        if a[0] == b[0]+b[3] or a[0]+a[3] == b[0]:
            dy = min(abs(a[1]+a[4] - b[1]), abs(b[1]+b[4] - a[1]))
            dz = min(abs(a[2]+a[5] - b[2]), abs(b[2]+b[5] - a[2]))

            if dy*dz==0:
                return False
            else:
                return float(dy*dz), float(d)

        elif a[1] == b[1]+b[4] or a[1]+a[4] == b[1]:
            dx = min(abs(a[0]+a[3] - b[0]), abs(b[0]+b[3] - a[0]))
            dz = min(abs(a[2]+a[5] - b[2]), abs(b[2]+b[5] - a[2]))

            if dx*dz==0:
                return False
            else:
                return float(dx*dz), float(d)

        elif a[2] == b[2]+b[5] or a[2]+a[5] == b[2]:
            dx = min(abs(a[0]+a[3] - b[0]), abs(b[0]+b[3] - a[0]))
            dy = min(abs(a[1]+a[4] - b[1]), abs(b[1]+b[4] - a[1]))

            if dx*dy==0:
                return False
            else:
                return float(dx*dy), float(d)

    else:
        return False

### Get all neighbours of a compartment c
### cmpts is the dictionary of all the compartments
def getNeighbours(cname, cdim, cmpts):
    nbrs = {}
    for nk,nv in cmpts.items():
        if not cname==nk:
            area = isNeighbour(np.array(cdim), np.array(nv))
            if area:
                nbrs.update({nk: area})
    return nbrs


### Get all vertices of a compartment
#@jit
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

    return np.array(verts)

### Draw each compartment of the bouton for visual verification
def plotCompartments(cmpts, bb):
    fig = plt.figure()
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
    plt.show()


### Make a list of initial index of each compartment
def initialIndex(cModels):
    cIdx = {}
    i = 0
    for cname, cm in cModels.items():
        #print cname, cModels[cname].nVar
        cIdx.update({cname: i})
        i += cm.nVar
    return cIdx

### Command line arguments
def commandArg(argv):
    for a in argv[1:]:
        a = a.split('=')
        cmdArg.update({a[0]: float(a[1])})


### get cModels containing the model equations
def getModels(cmpts, cHH, cVDCC, cPMCA, cPMCASensor):
    cModels = od()
    for cname, cdim in cmpts.items():
        if cname in cHH: model = mHH
        elif cname in cVDCC: model = mVDCC
        elif cname in cPMCA: model = mPMCA
        elif cname in cPMCASensor: model = mPMCASensor
        else: model = mCalbindin
        cModels.update({cname: mfb(model,
                                   name = cname,
                                   dim = cdim)
                      })
    return cModels

### A fancy progress bar
class FancyBar(IncrementalBar):
    t = 0
    suffix = Fore.CYAN + '%(done)0.2f msec' + Style.RESET_ALL
    def nextstep(self, t):
        self.t = t
        self.next()

    @property
    def done(self):
        return self.t


### The solution class
class solution:
    data = od()
    t = 0

    def __init__(self, cModels, cmpts, cmdArg, simName='trial/'):
        self.cModels = cModels
        self.cmdArg = cmdArg
        self.simName = simName

        for cname, cm in cModels.items():
            self.data.update({cname: cm.idx})
            cm.nbrs = getNeighbours(cm.name, cm.dim, cmpts)

    ### Putting all compartments together
    def dXdt(self, t, X):
        if t>self.t:

            self.bar.nextstep(1000*self.t)
            self.t += self.cmdArg['tf']/100.0

        #print 't =', self.t, 's'
        dX = []
        j=0
        for cm in self.cModels.values():
            ## All compartments have V value of (0,0,0)th compartment
            #cm.V = cModels[0].V

            ## Calcium Flux
            caFlux = 0
            for nbr, val in cm.nbrs.items():
                area, d = val
                diccCa = 1
                caFlux += diffCa*area*(X[self.cIdx[nbr]] - X[j])/d
                #if cm.name == '1-0-0':
                    #print t, cm.name, area, d, diffCa, (X[self.cIdx[nbr]] - X[j]), caFlux

            CaX = X[j:j+cm.nVar]
            CaX[0] += caFlux

            ## Collect dX values from each compartment
            dX += cm.dXdt(CaX, t)

            j += cm.nVar # increment counter to 1st element of next compartment

        return dX

    ### Solving the ODEs
    def solve(self):
        ## Make a list of initial index of each compartment
        self.cIdx = initialIndex(self.cModels)

        ## Simulation time
        ti, tf = 0, self.cmdArg['tf']
        tstep = self.cmdArg['tstep']
        self.tcp = self.cmdArg['tcp']

        ## initial values
        X0 = []
        for cm in self.cModels.values():
            X0 += cm.X0
        print 'Total number of equations:', Fore.GREEN, len(X0), Style.RESET_ALL, '\n'

        ## Solve ODE
        self.bar = FancyBar('Solving' + Fore.RED, max=101)

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
                for cname, cm in self.cModels.items():
                    header = 't\t\t'

                    for vn in self.data[cname]:
                        header += vn + '\t\t'
                        print vn

                    v = concatenate(([sol.t], sol.y[self.cIdx[cname]:self.cIdx[cname]+cm.nVar])).T
                    if temp:
                        savetxt(file, v, header=header, fmt='%.4e', delimiter='\t')
                    else:
                        savetxt(file, v, fmt='%.4e', delimiter='\t')

            ## Adding sol to solution till previous checkpoint
            if temp:
                t = sol.t
                y = sol.y
                temp = 0
            else:
                t = np.concatenate((t, sol.t))
                y = np.concatenate((y, sol.y), axis=1)

        self.bar.nextstep(1000*tf)
        self.bar.finish()

        if cmdArg['save']:
            file.close()
        ## Organise data in result.data dictionary
        for cname, idx in self.cIdx.items():
            for vname, v in self.data[cname].items():
                self.data[cname][vname] = y[idx+v]
        self.t = t
