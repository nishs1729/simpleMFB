from parameters import *

### Command line arguments
def commandArg(argv):
    for a in argv[1:]:
        a = a.split('=')
        cmdArg.update({a[0]: float(a[1])})


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
