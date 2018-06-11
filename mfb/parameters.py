from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from collections import OrderedDict as od
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
try:
    from colorama import Fore, Back, Style
except ImportError:
    print 'install <colorama> module:\n\tsudo pip install colorama\n'

from scipy import arange, linspace
import matplotlib.pyplot as plt
from scipy.integrate import *
import operator as op
from numba import *
import numpy as np
import sys, time
import warnings
import os

### Geometrical arrangement of all the compartments
### a:n   = a, a+1, a+2,...,n-1
### a:n:i = a, a+i, a+2*i,..., until (n-1)

### 2.4x2.4x1.4 um^3 with 200 nm as unit distance
modelDesc200 = '''[0:12, 0:12, 0]
                  [0:12:12, 0:12:12, 6]
                  [0:12:12,  0, 1:6:5]
                  [0:12:12, 11, 1:6:5]
                  [0,  1:11:10, 1:6:5]
                  [11, 1:11:10, 1:6:5]
                  [1:11:2, 1:11:2, 1:3:2]
                  [1:11:2, 1:11:2, 3:6:3]
                  unit=200
                  '''

### 2.4x2.4x1.4 um^3 with 100 nm as unit distance
modelDesc100 = '''[0:24, 0:24, 0]
                  [0:24:24, 0:24:24, 13]
                  [0:24:24,  0, 1:13:3]
                  [0:24:24, 23, 1:13:3]
                  [0,  1:23:22, 1:13:3]
                  [23, 1:23:22, 1:13:3]
                  [1:23:2, 1:23:2, 1]
                  [1:23:2, 1:23:2, 2:6:2]
                  [1:23:2, 1:23:2, 2:6:2]
                  [1,  1:23:22, 6:13:7]
                  [22, 1:23:22, 6:13:7]
                  [2:22:20,  1, 6:13:7]
                  [2:22:20, 22, 6:13:7]
                  [2:22:4, 2:22:4, 6:9:3]
                  [2:22:4, 2:22:4, 9:13:4]
                  unit=100
                  '''

### 2.4x2.4x1.4 um^3 with 50 nm as unit distance
modelDesc50 = '''[0:48, 0:48, 0]
                 [0:48:48, 0:48:48, 27]
                 [0:48:48,  0, 1:27:26]
                 [0:48:48, 47, 1:27:26]
                 [0,  1:47:46, 1:27:26]
                 [47, 1:47:46, 1:27:26]
                 [1:47:2, 1:47:2, 1:3:2]
                 [1:47:46, 46, 3:27:24]
                 [46, 1:46:45, 3:27:24]
                 [1:46:3, 1:46:3, 3:5:2]
                 [1:46:5, 1:46:5, 5:9:4]
                 [1:46:5, 1:46:5, 9:17:8]
                 [1:46:15, 1:46:15, 17:27:10]
                 unit=50
                 '''

### 4x2x1 um^3 with 0.1 um as unit distance
modelDesc421_100 = '''[0:2:2, 0:20:4, 0:2:3]
                      [38:40:2, 0:20:4, 0:2:3]
                      [0:40:5, 0:20:5, 7:10:3]
                      [0:40:4, 0:20:5, 5:7:2]
                      [0:40:2, 0:20:4, 3:5:2]
                      [2:38:2, 0:20:2, 2:3]
                      [2:38, 0:20, 0:2]
                      unit=100
                      '''

test1 = "[0,0,0]\n unit=200"
test8 = "[0:2,0:2,0:2]\n unit=200"


### Steady state initial values
initVal = {
    'Ca':        [100e-9],
    'HH':        [-64.0, 0.05, 0.6, 0.32],
    'PMCA':      [2.39e-6, 5.82e-7, 4.54e-10], # Total conc. = 2.98e-6 uM
    'caSensor':  [1.65e-6] + [0.0]*17,
    'VDCC':      [80e-6, 0.0, 0.0, 0.0, 0.0],
    'calbindin': [1.48e-05, 7.00e-06, 8.27e-07, 1.21e-05, 5.74e-06, 6.79e-07,
                  2.49e-06, 1.18e-06, 1.39e-07] # Total conc. = 45e-6 uM
}

### Command line arguments
cmdArg = od([
    ('tf',      100e-3),
    ('tstep',   1e-4),
    ('tcp',     50e-3),
    ('fig',     0),
    ('geo',     0),
    ('save',    0),
    ('dir',     'trial'),
    ('vfile',   0),
    ('rtol',    1e-4),
    ('atol',    1e-10),
    ('bar',     1),
    ('unit',    50) # smallest size unit in nm
])

try:
    from progress.bar import IncrementalBar
except ImportError:
    print 'install <progress> module for progress bar:'\
          + Fore.BLUE + '\n\tsudo pip install progress\n' + Fore.RESET
    cmdArg['bar'] = 0

### Parameters
## Calcium diffusion constant
diffCa = 2.2e-4 # m^2/s

## HH
C_m  = 1.0
g_Na = 120.0
g_K  = 36.0
g_L  = 0.3
E_Na = 50.0
E_K  = -77.0
E_L  = -54.387

## PMCA reaction rates
kPMCA01    = 1.5e8
kPMCA10    = 20
kPMCA12    = 100
kPMCA20    = 1e5
kPMCA0leak = 12.5

## Calbindin reaction rates
cbHon  = 0.55e7
cbHoff = 2.6
cbMon  = 4.35e7
cbMoff = 35.8

## VDCC
a10, a20, a30, a40 = 4040, 6700, 4390, 17330 # /sec
b10, b20, b30, b40 = 2880, 6300, 8160, 1840  # /sec
V1,  V2,  V3,  V4  = 49.14, 42.08, 55.31, 26.55 # mV

## Calcium Sensor
sf = 0.612e8    # /M s
sb = 2.32e3     # /s
af = 3.82e6     # /s
ab = 13         # /s
b  = 0.25       # /s
rsy = 0
ras = 0
rsp = 0

### Compartment list for specific model type
cm = {
    'cHH': ['0-0-0'],
    'cVDCC': ['1-0-0', '1-1-0', '1-0-1', '1-1-1', '0-0-1', '0-1-0', '0-1-1'],
    'cPMCA': [],
    'cPMCASensor': []
}

### All available model components
models = {
    'Ca': [],
    'HH': [],
    'PMCA': [],
    'VDCC': [],
    'caSensor': [],
    'calbindin': []
}

### compartment with PMCA and VDCC models
mVDCC = {
    'Ca': [],
    'PMCA': [],
    #'VDCC': [],
    #'calbindin': []
}

### compartment with PMCA model
mPMCA = {
    'Ca': [],
    'PMCA': [],
    'calbindin': []
}

### compartment with PMCA and caSensor models
mPMCASensor = {
    'Ca': [],
    'PMCA': [],
    'caSensor': [],
    'calbindin': []
}

### compartment with HH model (0-0-0)
mHH = {
    'Ca': [],
    'HH': [],
    'PMCA': [],
    #'calbindin': []
}

### compartment with Calbindin model
mCalbindin = {
    'Ca': [],
    'PMCA': []
    #'calbindin': []
}
#set_printoptions(precision=4)
