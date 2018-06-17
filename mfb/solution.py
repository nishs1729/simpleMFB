from parameters import *
from MFBFunctions import *
try:
    from misc import FancyBar
except ImportError:
    cmdArg['bar'] = 0

from misc import getV

### The solution class
class solution:
    data = od()
    t = 0

    def __init__(self, cModels, cmpts, simName='trial/'):
        self.cModels = cModels
        self.simName = simName

        if cmdArg['vfile']:
            self.vfile = getV('v.txt')

        for cname, cm in cModels.items():
            self.data.update({cname: cm.idx})
            cm.nbrs = getNeighbours(cm.name, cm.dim, cmpts)

    ### Putting all compartments together
    def dXdt(self, X, t):
        if t>self.t:
            if cmdArg['bar']:
                self.bar.nextstep(1000*self.t, time.time()-self.timei)
                self.t += cmdArg['tf']/100.0
            else:
                self.t += cmdArg['tf']/20
                print 't:', self.t*1000, 'msec'

        dX = []
        j=0
        for cm in self.cModels.values():
            ## All compartments have V value of (0,0,0)th compartment
            #cm.V = cModels['0-0-0'].V
            if cmdArg['vfile']:
                cm.V = self.vfile(t)

            ## Calcium Flux
            caFlux = 0
            for nbr, val in cm.nbrs.items():
                area, d = val
                caFlux += diffCa*area*(X[self.cIdx[nbr]] - X[j])/d/cm.vol*1e18

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
        self.timei = time.time()

        ## Simulation time
        ti, tf = 0, cmdArg['tf']
        tstep = cmdArg['tstep']
        self.tcp = cmdArg['tcp']

        ## initial values
        X0 = []
        for cm in self.cModels.values():
            X0 += cm.X0
        print 'Total equations:' + Fore.GREEN, len(X0), Style.RESET_ALL, '\n'

        ## Initialise Progressbar
        if cmdArg['bar']:
            self.bar = FancyBar('Solving' + Fore.RED, max=101)

        tinterval = []
        aa = arange(ti, tf, self.tcp)
        for i in range(len(aa[:-1])):
            tinterval += [[aa[i], aa[i+1]]]
        tinterval += [[aa[-1], tf]]

        ## Saving data till current checkpoint in files
        if cmdArg['save']:
            dir = 'data/'+cmdArg['dir']+'/'
            if not os.path.exists(dir):
                os.makedirs(dir)
            file = {}
            for cname, cm in self.cModels.items():
                fname = dir+cname+'.txt'
                file.update({cname: open(fname, 'w')})

            ## Save the X0 values at each checkpoint
            fname = dir+'.chkptX0'
            chkptX0 = open(fname, 'w')


        temp = 1
        for ti, tf in tinterval:
            t_eval = linspace(ti, tf, round((tf - ti)/tstep) + 1)[:-1]

            sol = odeint(self.dXdt, X0, t_eval)

            ## Last values of sol as X0
            X0 = sol[-1,:]

            ## Saving data till current checkpoint in files
            if cmdArg['save']:
                for cname, cm in self.cModels.items():
                    header = 't\t\t'

                    for vn in self.data[cname]:
                        header += vn + '\t\t'

                    v = np.concatenate(([t_eval], sol.T[self.cIdx[cname]:self.cIdx[cname]+cm.nVar])).T
                    if temp:
                        np.savetxt(file[cname], v, header=header, fmt='%.4e', delimiter='\t')
                    else:
                        np.savetxt(file[cname], v, fmt='%.4e', delimiter='\t')

                ## Save X0 at each checkpoint in file '.chkptX0'
                np.savetxt(chkptX0, [[t_eval[-1]] + X0.tolist()], fmt='%.4e', delimiter='\t')

            ## Adding sol to solution till previous checkpoint
            if temp:
                t = t_eval
                y = sol.T
                temp = 0
            else:
                t = np.concatenate((t, t_eval))
                y = np.concatenate((y, sol.T), axis=1)

        if cmdArg['bar']:
            self.bar.nextstep(1000*tf, time.time()-self.timei)
            self.bar.finish()

        if cmdArg['save']:
            for cname in self.cModels.keys():
                file[cname].close()

        ## Organise data in result.data dictionary
        for cname, idx in self.cIdx.items():
            for vname, v in self.data[cname].items():
                self.data[cname][vname] = y[idx+v]
        self.t = t
