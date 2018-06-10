from parameters import *

### Command line arguments
def commandArg(argv):
    for a in argv[1:]:
        ## See default arguments
        if a == "--default" or a == "-d":
            print '\nDefault command arguments:\n'
            for k,v in cmdArg.items():
                print k, ':\t', v
            exit()

        ## Update command arguments
        a = a.split('=')
        if a[0] in ('fig', 'geo', 'save', 'vfile'): ##for int
            cmdArg.update({a[0]: int(a[1])})
        elif a[0] in ('tf', 'tstep', 'tcp', 'rtol', 'atol'): ## for float
            cmdArg.update({a[0]: float(a[1])})
        elif a[0] in ('dir'): ## for directory
            cmdArg.update({a[0]: a[1]+'/'})
        else: ## for string
            cmdArg.update({a[0]: a[1]})

### A fancy progress bar
class FancyBar(IncrementalBar):
    t_sim = 0
    t_real = 0
    suffix = Fore.CYAN + '%(simTime)0.2f msec' + Fore.RED + ' [Real Time:'\
             + '%(realTime)d sec]' + Style.RESET_ALL
    def nextstep(self, t_sim, t_real):
        self.t_sim = t_sim
        self.t_real = t_real
        self.next()

    @property
    def simTime(self):
        return self.t_sim

    @property
    def realTime(self):
        return self.t_real


### decorator function to time functions
def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        tf = time.time()
        print 'Func: ' + method.__name__ + ' run in ' + str(tf-ts) + 's'
        return result
    return timed

### Get timeseries voltage from file and return interpolated function
def getV(fname='v.txt'):
    data = np.genfromtxt(fname, unpack=True, usecols=(0,1))
    return interp1d(data[0], data[1], kind='cubic')
