from parameters import *

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
def b1(V):    return b10*np.exp(-V/V1)

def a2(V):    return a20*np.exp(V/V2)
def b2(V):    return b20*np.exp(-V/V2)

def a3(V):    return a30*np.exp(V/V3)
def b3(V):    return b30*np.exp(-V/V3)

def a4(V):    return a40*np.exp(V/V4)
def b4(V):    return b40*np.exp(-V/V4)
