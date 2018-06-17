import numpy as np

### HH model
C_m  = 1.0
g_Na = 120.0
g_K  = 36.0
g_L  = 0.3
E_Na = 50.0
E_K  = -77.0
E_L  = -54.387

### External Current for HH model
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


### pqVDCC
pq_a10, pq_a20, pq_a30, pq_a40, pq_a = 5890, 9210, 5200, 1823180, 247710 # /sec
pq_b10, pq_b20, pq_b30, pq_b40, pq_b = 14990, 6630, 132800, 248580, 8280 # /sec
pq_V1,  pq_V2,  pq_V3,  pq_V4        = 62.61, 33.92, 135.08, 20.86 # mV

### pqVDCC gating variables
def pq_a1(V):    return pq_a10*np.exp( V/pq_V1)
def pq_b1(V):    return pq_b10*np.exp(-V/pq_V1)
def pq_a2(V):    return pq_a20*np.exp( V/pq_V2)
def pq_b2(V):    return pq_b20*np.exp(-V/pq_V2)
def pq_a3(V):    return pq_a30*np.exp( V/pq_V3)
def pq_b3(V):    return pq_b30*np.exp(-V/pq_V3)
def pq_a4(V):    return pq_a40*np.exp( V/pq_V4)
def pq_b4(V):    return pq_b40*np.exp(-V/pq_V4)


### rVDCC
r_a10, r_a20, r_a30, r_a40, r_a =  9911360, 4880, 4000, 256410, 228830 # /sec
r_b10, r_b20, r_b30, r_b40, r_b =  620, 21910, 51300, 116970, 1780 # /sec
r_V1,  r_V2,  r_V3,  r_V4       =  67.75, 50.94, 173.29, 16.92 # mV

### rVDCC gating variables
def r_a1(V):    return r_a10*np.exp( V/r_V1)
def r_b1(V):    return r_b10*np.exp(-V/r_V1)
def r_a2(V):    return r_a20*np.exp( V/r_V2)
def r_b2(V):    return r_b20*np.exp(-V/r_V2)
def r_a3(V):    return r_a30*np.exp( V/r_V3)
def r_b3(V):    return r_b30*np.exp(-V/r_V3)
def r_a4(V):    return r_a40*np.exp( V/r_V4)
def r_b4(V):    return r_b40*np.exp(-V/r_V4)


### nVDCC
n_a10, n_a20, n_a30, n_a40, n_a = 4290, 5240, 4980, 772630, 615010 # /sec
n_b10, n_b20, n_b30, n_b40, n_b = 5230, 6630, 73890, 692180, 7680 # /sec
n_V1,  n_V2,  n_V3,  n_V4       = 68.75, 39.53, 281.62, 18.46 # mV

### nVDCC gating variables
def n_a1(V):    return n_a10*np.exp( V/n_V1)
def n_b1(V):    return n_b10*np.exp(-V/n_V1)
def n_a2(V):    return n_a20*np.exp( V/n_V2)
def n_b2(V):    return n_b20*np.exp(-V/n_V2)
def n_a3(V):    return n_a30*np.exp( V/n_V3)
def n_b3(V):    return n_b30*np.exp(-V/n_V3)
def n_a4(V):    return n_a40*np.exp( V/n_V4)
def n_b4(V):    return n_b40*np.exp(-V/n_V4)


### pqVDCC from CA3 MCell model
#pq_a10, pq_a20, pq_a30, pq_a40 = 4040, 6700, 4390, 17330 # /sec
#pq_b10, pq_b20, pq_b30, pq_b40 = 2880, 6300, 8160, 1840  # /sec
#pq_V1,  pq_V2,  pq_V3,  pq_V4  = 49.14, 42.08, 55.31, 26.55 # mV
"""
### pqVDCC gating variables from CA3 MCell model
def a1(V):    return pq_a10*np.exp( V/pq_V1)
def b1(V):    return pq_b10*np.exp(-V/pq_V1)
def a2(V):    return pq_a20*np.exp( V/pq_V2)
def b2(V):    return pq_b20*np.exp(-V/pq_V2)
def a3(V):    return pq_a30*np.exp( V/pq_V3)
def b3(V):    return pq_b30*np.exp(-V/pq_V3)
def a4(V):    return pq_a40*np.exp( V/pq_V4)
def b4(V):    return pq_b40*np.exp(-V/pq_V4)

### pqVDCC equation from CA3 MCell model
if 'pqVDCC' in self.models:
    dCa += 19.3*self.V*(0.3993 - np.exp(-self.V/80.36))/(1 - np.exp(self.V/80.36))*pqVDCC_O
    dpqVDCC_C0 = + pq_b1(self.V)*pqVDCC_C1 - pq_a1(self.V)*pqVDCC_C0
    dpqVDCC_C1 = + pq_a1(self.V)*pqVDCC_C0 + pq_b2(self.V)*pqVDCC_C2 \
                 -(pq_b1(self.V) + pq_a2(self.V))*pqVDCC_C1
    dpqVDCC_C2 = + pq_a2(self.V)*pqVDCC_C1 + pq_b3(self.V)*pqVDCC_C3 \
                 -(pq_b2(self.V) + pq_a3(self.V))*pqVDCC_C2
    dpqVDCC_C3 = + pq_a3(self.V)*pqVDCC_C2 + pq_b4(self.V)*pqVDCC_O \
                 -(pq_b3(self.V) + pq_a4(self.V))*pqVDCC_C3
    dpqVDCC_O  = + pq_a4(self.V)*pqVDCC_C3 - pq_b4(self.V)*pqVDCC_O

    self.dX += [dpqVDCC_C0, dpqVDCC_C1, dpqVDCC_C2, dpqVDCC_C3, dpqVDCC_O]
    self.dX[0] += dCa
"""
