from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from collections import OrderedDict as od
from mpl_toolkits.mplot3d import Axes3D
from progress.bar import IncrementalBar
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
    'tstep': 1e-3,
    'fig'  : 0,
    'rtol' : 1e-4,
    'atol' : 1e-10,
    'simName': 'trial/'
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

## for division by zero cases
delta = 1e-4

### External Current
def I_inj(t):  return 10#*(t>0.005) - 10*(t>0.015) + 35*(t>0.3) - 35*(t>0.4)

### Channel gating variables (ms)
def alpha_m(V):   return 0.1*(V+40.0)/(1.0 - np.exp(-(V+40.0) / 10.0))
def beta_m(V):    return 4.0*np.exp(-(V+65.0) / 18.0)
def alpha_h(V):   return 0.07*np.exp(-(V+65.0) / 20.0)
def beta_h(V):    return 1.0/(1.0 + np.exp(-(V+35.0) / 10.0))
def alpha_n(V):   return 0.01*(V+55.0)/(1.0 - np.exp(-(V+55.0) / 10.0))
def beta_n(V):    return 0.125*np.exp(-(V+65.0) / 80.0)

### Membrane current (in uA/cm^2)
def I_Na(V, m, h):  return g_Na * m**3 * h * (V - E_Na)
def I_K(V, n):      return g_K  * n**4 * (V - E_K)
def I_L(V):         return g_L * (V - E_L)

### VDCC gating variables
def a1(V):    return a10*np.exp(V/V1)
def b1(V):    return b10*np.exp(-V/V1)#!/usr/bin/env python

def a2(V):    return a20*np.exp(V/V2)
def b2(V):    return b20*np.exp(-V/V2)

def a3(V):    return a30*np.exp(V/V3)
def b3(V):    return b30*np.exp(-V/V3)

def a4(V):    return a40*np.exp(V/V4)
def b4(V):    return b40*np.exp(-V/V4)

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
