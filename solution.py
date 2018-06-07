from MFBfunctions import *
from parameters import *

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
        if t>=self.t:
            #print 't =', self.t, 's'
            self.bar.nextstep(self.t)
            self.t += self.tcp/5

        dX = []
        j=0
        for cm in self.cModels.values():
            ## All compartments have V value of (0,0,0)th compartment
            #cm.V = cModels[0].V

            ## Calcium Flux
            caFlux = 0
            for nbr, val in cm.nbrs.items():
                area, d = val
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
        ti, tf = 0.0, self.cmdArg['tf']
        tstep = self.cmdArg['tstep']
        self.tcp = self.cmdArg['tcp']

        ## initial values
        X0 = []
        for cm in self.cModels.values():
            X0 += cm.X0
        print 'Total number of equations:', len(X0)

        ## Solve ODE
        print "solving the ODE...", tf/self.tcp*5
        self.bar = FancyBar('Solving', max=int(tf/self.tcp*5))

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

        self.bar.finish()

        if cmdArg['save']:
            file.close()
        ## Organise data in result.data dictionary
        for cname, idx in self.cIdx.items():
            for vname, v in self.data[cname].items():
                self.data[cname][vname] = y[idx+v]
        self.t = t
