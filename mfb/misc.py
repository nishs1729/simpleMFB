from parameters import *

### Command line arguments
def commandArg(argv):
    for a in argv[1:]:
        a = a.split('=')
        cmdArg.update({a[0]: float(a[1])})


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
