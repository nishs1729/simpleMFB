from parameters import *
from modelEquations import equations

### Takes in location and size of each compartment
### and returns a dictionary of compartment name and
### corresponding location ad size array
def compartments(modelDesc):
    desc = np.array([a.strip() for a in modelDesc.strip().split("\n")])
    cmpts = od()
    cSurf = []
    for a in desc:
        if 'unit' in a: # neglect the unit size
            continue

        a = a.strip('[').strip(']').split(',')
        #print 'a:', a

        arr, temp = [], []
        for i in a[:3]:
            if ':' in i:
                t = np.array([int(a) for a in i.split(':')])
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

### Return the bounding box for the compartments
def boundingBox(desc):
    desc = np.array([a.strip() for a in desc.strip().split("\n")])

    size, size1 = [0, 0, 0], [0, 0, 0]
    flag = [0, 0, 0]
    for a in desc:
        if 'unit' in a: # update the unit size
            a = int(a.split('=')[1])
            cmdArg['unit'] = a
            continue

        a = a.strip('[').strip(']').split(',')

        for i in range(len(size)):
            size[i] = max([int(s) for s in a[i].split(':')] + size[i:i+1])
            if ':' not in a[i]:
                size1[i] = max(size1[i], int(a[i]))
            if size[i] > size1[i]:
                flag[i] = 0
            else:
                flag[i] = 1

    bbmax = [str(a) for a in map(lambda x, y: x+y, size, flag)]
    bBox = "[0:" + bbmax[0] + ", 0:" + bbmax[1] + ", 0:" + bbmax[2] + "]"
    return bBox

### Check if no compartment have overlapping volumes and
### there are no gaps in the model.
### comp is the bounding box of all compartments broken into
### smallest cubic unit. For eg. "[0:3,0:3,0:3]".
### cmpts is the dictionary of all the compartments.
def checkGeometry(bBox, cmpts):
    sCubes = compartments(bBox)

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
        print Fore.RED + 'Error in compartment geometry!'
        print 'gaps:', gap, len(gap)
        print 'overlaps:', overlap, Fore.RESET
        return False


### Check if two compartments a amd b are neighbours
### Returns False if they are not neighbours
### Returns shared surface area if they are neighbours
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
        u = cmdArg['unit']

        ## If a and b have an overlapping surface, return area and d
        if a[0] == b[0]+b[3] or a[0]+a[3] == b[0]:
            dy = min(abs(a[1]+a[4] - b[1]), abs(b[1]+b[4] - a[1]))
            dz = min(abs(a[2]+a[5] - b[2]), abs(b[2]+b[5] - a[2]))

            if dy*dz==0:
                return False
            else:
                return float(dy*dz*u*u), float(d*u)

        elif a[1] == b[1]+b[4] or a[1]+a[4] == b[1]:
            dx = min(abs(a[0]+a[3] - b[0]), abs(b[0]+b[3] - a[0]))
            dz = min(abs(a[2]+a[5] - b[2]), abs(b[2]+b[5] - a[2]))

            if dx*dz==0:
                return False
            else:
                return float(dx*dz*u*u), float(d*u)

        elif a[2] == b[2]+b[5] or a[2]+a[5] == b[2]:
            dx = min(abs(a[0]+a[3] - b[0]), abs(b[0]+b[3] - a[0]))
            dy = min(abs(a[1]+a[4] - b[1]), abs(b[1]+b[4] - a[1]))

            if dx*dy==0:
                return False
            else:
                return float(dx*dy*u*u), float(d*u)

    else:
        return False

### Get all neighbours of a compartment c
### cmpts is the dictionary of all the compartments
def getNeighbours(cname, cdim, cmpts):
    nbrs = {}
    for nk,nv in cmpts.items():
        if not cname==nk:
            area_d = isNeighbour(np.array(cdim), np.array(nv))
            if area_d:
                nbrs.update({nk: area_d}) # area in nm^2 | d in nm
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

    return np.array(verts)

### Draw each compartment of the bouton for visual verification
def plotCompartments(cmpts, bBox, cm):
    bBox = bBox.strip('[').strip(']').split(',')
    bBox = [float(a) for a in [bb.split(':')[1] for bb in bBox]]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # set x,y,z limits
    ax.set_xlim(0, bBox[0])
    ax.set_ylim(0, bBox[1])
    ax.set_zlim(0, bBox[2])

    ax.axis('off')

    # plot sides
    for c in cmpts.values():
        verts = getVertices(c)
        collection = Poly3DCollection(verts, linewidths=0.5,
                                    edgecolors='g', alpha=.2)
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


### get cModels containing the model equations
def getModels(cmpts, cm, initC):
    model = {}
    ## All compartments have Ca
    for m in cmpts.keys():
        if m in initC.keys():
            if 'Ca' in initC[m].keys():
                model.update({m: {'Ca': initC[m]['Ca']}})
        else:
            model.update({m: {'Ca': []}})

    ## Fill compartments with models as required
    for m, cnames in cm.items():
        for cname in cnames:
            try:
                if m in initC[cname].keys():
                    model[cname].update({m: initC[cname][m]})
                else:
                    model[cname].update({m: []})
            except:
                model[cname].update({m: []})

    ## Generate models for compartments
    cModels = od()
    for cname, cdim in cmpts.items():
        cModels.update({cname: equations(model[cname],
                               name = cname,
                               dim = cdim)})

    return cModels


### Return names of components containing requested points
### length, l : nm
### center    : nm
def hexPoints(num, l, bBox, center=[0, 0]):
    shift = [np.array(a) for a in [[l,0], [0,l], [-l,0], [0,-l]]]
    points = [np.array(center)]
    sh, i, n = 0, 1, 1
    while True:
        br = False
        for _ in range(2):
            for __ in range(n):
                i += 1
                if i>num:
                    br = True
                    break
                points.append(points[-1] + shift[sh])
            sh = (sh+1)%4
            if br:
                break
        n += 1
        if br:
            break

    ### Get compartments
    u = cmdArg['unit']
    xy = np.array([aa.split(':') for aa in bBox.strip('[]').split(',')[:2]]).flatten()
    xmin, xmax, ymin, ymax = [int(a) for a in xy]
    midx, midy = int((xmax-xmin)/2.0), int((ymax-ymin)/2.0)

    cname = []
    for px, py in points:
        #print px + midx*u, py + midy*u, '\t', int(midx + px/u), int(midy + py/u)
        x, y = str(int(midx + px/u)), str(int(midy + py/u))
        cname.append(x + '-' + y + '-0')

    return cname

### Return names of surface compartments
def cSurf(cmpts):
    pa
