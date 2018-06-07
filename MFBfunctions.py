from parameters import *

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


### A fancy progress bar
class FancyBar(ChargingBar):
    t = 0
    suffix = '%(done).2f sec'
    def nextt(self, t):
        self.t = t
        self.next()

    @property
    def done(self):
        return self.t
