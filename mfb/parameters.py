from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from collections import OrderedDict as od
from mpl_toolkits.mplot3d import Axes3D
from progress.bar import IncrementalBar
from scipy.interpolate import interp1d
from colorama import Fore, Back, Style
from scipy import arange, linspace
import matplotlib.pyplot as plt
from scipy.integrate import *
import operator as op
from numba import *
import numpy as np
import sys, time
import warnings
import os

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
cmdArg = {
    'tcp'  : 50e-3,
    'geo'  : 0,
    'save' : 0,
    'tf'   : 100e-3,
    'tstep': 1e-4,
    'fig'  : 0,
    'rtol' : 1e-4,
    'atol' : 1e-10,
    'simName': 'trial/',
    'vfile': 0
}

### Parameters
## Calcium diffusion constant
diffCa = 220

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
