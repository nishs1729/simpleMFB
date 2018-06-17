from collections import OrderedDict as od
from modelFunctions import *
from parameters import *

class equations:
    def __init__(self, models, name, dim):
        self.models = models
        self.name = name
        self.dim  = dim
        self.vol  = reduce(lambda x, y: x*y, dim[3:])*cmdArg['unit']**3
        #print nbrs, dim, self.vol

        # Indexing of the compartment variables
        i = 0
        self.idx = od()
        self.X0 = []
        if 'Ca' in models:
            if models['Ca'] == []: self.X0 += initVal['Ca']
            else: self.X0 += models['Ca']
            self.idx.update(od([('Ca', i)]))
            i += len(initVal['Ca'])

        if 'HH' in models:
            if models['HH'] == []: self.X0 += initVal['HH']
            else: self.X0 += models['HH']
            self.idx.update(od(zip(['V', 'm', 'h', 'n'], range(i,i+4))))
            i += len(initVal['HH'])
            self.V = 0

        if 'PMCA' in models:
            if models['PMCA'] == []: self.X0 += initVal['PMCA']
            else: self.X0 += models['PMCA']
            self.idx.update(od(zip(['PMCA0', 'PMCA1', 'PMCA2'], range(i,i+3))))
            i += len(initVal['PMCA'])

        if 'pqVDCC' in models:
            if models['pqVDCC'] == []: self.X0 += initVal['pqVDCC']
            else: self.X0 += models['pqVDCC']
            self.idx.update(od(zip(['pqVDCC_C0', 'pqVDCC_C1', 'pqVDCC_C2', 'pqVDCC_C3', 'pqVDCC_C4', 'pqVDCC_O'],
                                    range(i,i+6))))
            i += len(initVal['pqVDCC'])
            self.V = 0

        if 'rVDCC' in models:
            if models['rVDCC'] == []: self.X0 += initVal['rVDCC']
            else: self.X0 += models['rVDCC']
            self.idx.update(od(zip(['rVDCC_C0', 'rVDCC_C1', 'rVDCC_C2', 'rVDCC_C3', 'rVDCC_C4', 'rVDCC_O'],
                                    range(i,i+6))))
            i += len(initVal['rVDCC'])
            self.V = 0

        if 'nVDCC' in models:
            if models['nVDCC'] == []: self.X0 += initVal['nVDCC']
            else: self.X0 += models['nVDCC']
            self.idx.update(od(zip(['nVDCC_C0', 'nVDCC_C1', 'nVDCC_C2', 'nVDCC_C3', 'nVDCC_C4', 'nVDCC_O'],
                                    range(i,i+6))))
            i += len(initVal['nVDCC'])
            self.V = 0

        if 'calbindin' in models:
            if models['calbindin'] == []: self.X0 += initVal['calbindin']
            else: self.X0 += models['calbindin']
            self.idx.update(
                        od([('cbH0M0', i)  , ('cbH0M1', i+1), ('cbH0M2', i+2),
                            ('cbH1M0', i+3), ('cbH1M1', i+4), ('cbH1M2', i+5),
                            ('cbH2M0', i+6), ('cbH2M1', i+7), ('cbH2M2', i+8)]))
            i += len(initVal['calbindin'])

        if 'AZ' in models:
            if models['AZ'] == []: self.X0 += initVal['AZ']
            else: self.X0 += models['AZ']
            self.idx.update(
                    od([('AZ00', i+0 ), ('AZ10', i+1 ), ('AZ20', i+2 ),
                        ('AZ30', i+3 ), ('AZ40', i+4 ), ('AZ50', i+5 ),
                        ('AZ01', i+6 ), ('AZ11', i+7 ), ('AZ21', i+8 ),
                        ('AZ31', i+9 ), ('AZ41', i+10), ('AZ51', i+11),
                        ('AZ02', i+12), ('AZ12', i+13), ('AZ22', i+14),
                        ('AZ32', i+15), ('AZ42', i+16), ('AZ52', i+17)]))
            i += len(initVal['AZ'])

        self.nVar = i
        #rint 'X0:', self.X0, '\n', idx, '\n', self.nVar, '\n\n'


    def dXdt(self, X, t):
        self.dX = []
        i=0
        if 'Ca' in self.models:
            j = len(initVal['Ca'])
            Ca = X[i]
            i += j

        if 'HH' in self.models:
            j = len(initVal['HH'])
            self.V, m, h, n = X[i:i+j]  # V is a property of object as it is to be provided externally
            i += j

        if 'PMCA' in self.models:
            j = len(initVal['PMCA'])
            PMCA0, PMCA1, PMCA2 = X[i:i+j]
            i += j

        if 'calbindin' in self.models:
            j = len(initVal['calbindin'])
            cbH0M0, cbH0M1, cbH0M2, cbH1M0,\
            cbH1M1, cbH1M2, cbH2M0, cbH2M1, cbH2M2 = X[i:i+j]
            i += j

        if 'pqVDCC' in self.models:
            j = len(initVal['pqVDCC'])
            pqVDCC_C0, pqVDCC_C1, pqVDCC_C2, pqVDCC_C3, pqVDCC_C4, pqVDCC_O = X[i:i+j]
            i += j

        if 'rVDCC' in self.models:
            j = len(initVal['rVDCC'])
            rVDCC_C0, rVDCC_C1, rVDCC_C2, rVDCC_C3, rVDCC_C4, rVDCC_O = X[i:i+j]
            i += j

        if 'nVDCC' in self.models:
            j = len(initVal['nVDCC'])
            nVDCC_C0, nVDCC_C1, nVDCC_C2, nVDCC_C3, nVDCC_C4, nVDCC_O = X[i:i+j]
            i += j

        if 'AZ' in self.models:
            j = len(initVal['AZ'])
            AZ00, AZ10, AZ20, AZ30, AZ40, AZ50, AZ01, AZ11, AZ21, \
            AZ31, AZ41, AZ51, AZ02, AZ12, AZ22, AZ32, AZ42, AZ52 = X[i:i+j]
            i += j


        #################### Equations #####################
        if 'Ca' in self.models:
            dCa = 0
            self.dX += [dCa]

        ### HH model
        if 'HH' in self.models:
            dV = 1000*(I_inj(t) - I_Na(self.V, m, h) - I_K(self.V, n) - I_L(self.V)) / C_m
            dm = 1000*(alpha_m(self.V)*(1.0-m) - beta_m(self.V)*m)
            dh = 1000*(alpha_h(self.V)*(1.0-h) - beta_h(self.V)*h)
            dn = 1000*(alpha_n(self.V)*(1.0-n) - beta_n(self.V)*n)

            self.dX += [dV, dm, dh, dn]

        ### PMCA
        if 'PMCA' in self.models:
            dCa   += + (kPMCA0leak - kPMCA01*Ca)*PMCA0 + kPMCA10*PMCA1
            dPMCA0 = + kPMCA20*PMCA2 + kPMCA10*PMCA1 - kPMCA01*Ca*PMCA0
            dPMCA1 = - (kPMCA12 + kPMCA10)*PMCA1 + kPMCA01*Ca*PMCA0
            dPMCA2 = - kPMCA20*PMCA2 + kPMCA12*PMCA1

            self.dX += [dPMCA0, dPMCA1, dPMCA2]
            self.dX[0] += dCa

        ### pqVDCC
        if 'pqVDCC' in self.models:
            '''dCa += 0 #1*pqVDCC_O
            dpqVDCC_C0 = + pq_b1(self.V)*pqVDCC_C1 - pq_a1(self.V)*pqVDCC_C0
            dpqVDCC_C1 = + pq_a1(self.V)*pqVDCC_C0 + pq_b2(self.V)*pqVDCC_C2 \
                         -(pq_b1(self.V) + pq_a2(self.V))*pqVDCC_C1
            dpqVDCC_C2 = + pq_a2(self.V)*pqVDCC_C1 + pq_b3(self.V)*pqVDCC_C3 \
                         -(pq_b2(self.V) + pq_a3(self.V))*pqVDCC_C2
            dpqVDCC_C3 = + pq_a3(self.V)*pqVDCC_C2 + pq_b4(self.V)*pqVDCC_C4 \
                         -(pq_b3(self.V) + pq_a4(self.V))*pqVDCC_C3
            dpqVDCC_C4 = + pq_a4(self.V)*pqVDCC_C3 + pq_b*pqVDCC_O \
                         -(pq_b4(self.V) + pq_a)*pqVDCC_C4
            dpqVDCC_O  = + pq_a*pqVDCC_C4 - pq_b*pqVDCC_O

            self.dX += [dpqVDCC_C0, dpqVDCC_C1, dpqVDCC_C2, dpqVDCC_C3, dpqVDCC_C4, dpqVDCC_O]
            self.dX[0] += dCa
            '''
            ### pqVDCC equation from CA3 MCell model
            dCa += 19.3*self.V*(0.3993 - np.exp(-self.V/80.36))/(1 - np.exp(self.V/80.36))*pqVDCC_O
            dpqVDCC_C0 = + pq_b1(self.V)*pqVDCC_C1 - pq_a1(self.V)*pqVDCC_C0
            dpqVDCC_C1 = + pq_a1(self.V)*pqVDCC_C0 + pq_b2(self.V)*pqVDCC_C2 \
                         -(pq_b1(self.V) + pq_a2(self.V))*pqVDCC_C1
            dpqVDCC_C2 = + pq_a2(self.V)*pqVDCC_C1 + pq_b3(self.V)*pqVDCC_C3 \
                         -(pq_b2(self.V) + pq_a3(self.V))*pqVDCC_C2
            dpqVDCC_C3 = + pq_a3(self.V)*pqVDCC_C2 + pq_b4(self.V)*pqVDCC_O \
                         -(pq_b3(self.V) + pq_a4(self.V))*pqVDCC_C3
            dpqVDCC_C4 = 0
            dpqVDCC_O  = + pq_a4(self.V)*pqVDCC_C3 - pq_b4(self.V)*pqVDCC_O

            self.dX += [dpqVDCC_C0, dpqVDCC_C1, dpqVDCC_C2, dpqVDCC_C3, dpqVDCC_C4, dpqVDCC_O]
            self.dX[0] += dCa
            #'''

        ### rVDCC
        if 'rVDCC' in self.models:
            dCa += 1*rVDCC_O
            drVDCC_C0 = + r_b1(self.V)*rVDCC_C1 - r_a1(self.V)*rVDCC_C0
            drVDCC_C1 = + r_a1(self.V)*rVDCC_C0 + r_b2(self.V)*rVDCC_C2 \
                        -(r_b1(self.V) + r_a2(self.V))*rVDCC_C1
            drVDCC_C2 = + r_a2(self.V)*rVDCC_C1 + r_b3(self.V)*rVDCC_C3 \
                        -(r_b2(self.V) + r_a3(self.V))*rVDCC_C2
            drVDCC_C3 = + r_a3(self.V)*rVDCC_C2 + r_b4(self.V)*rVDCC_C3 \
                        -(r_b3(self.V) + r_a4(self.V))*rVDCC_C3
            drVDCC_C4 = + r_a4(self.V)*rVDCC_C3 + r_b*rVDCC_O \
                        -(r_b4(self.V) + r_a)*rVDCC_C4
            drVDCC_O  = + r_a*rVDCC_C4 - r_b*rVDCC_O

            self.dX += [drVDCC_C0, drVDCC_C1, drVDCC_C2, drVDCC_C3, drVDCC_C4, drVDCC_O]
            self.dX[0] += dCa

        ### nVDCC
        if 'nVDCC' in self.models:
            dCa += 1000*nVDCC_O
            dnVDCC_C0 = + r_b1(self.V)*nVDCC_C1 - r_a1(self.V)*nVDCC_C0
            dnVDCC_C1 = + r_a1(self.V)*nVDCC_C0 + r_b2(self.V)*nVDCC_C2 \
                         -(r_b1(self.V) + r_a2(self.V))*nVDCC_C1
            dnVDCC_C2 = + r_a2(self.V)*nVDCC_C1 + r_b3(self.V)*nVDCC_C3 \
                         -(r_b2(self.V) + r_a3(self.V))*nVDCC_C2
            dnVDCC_C3 = + r_a3(self.V)*nVDCC_C2 + r_b4(self.V)*nVDCC_C4 \
                         -(r_b3(self.V) + r_a4(self.V))*nVDCC_C3
            dnVDCC_C4 = + n_a4(self.V)*nVDCC_C3 + n_b*nVDCC_O \
                         -(n_b4(self.V) + n_a)*nVDCC_C4
            dnVDCC_O  = + n_a*nVDCC_C4 - n_b*nVDCC_O

            self.dX += [dnVDCC_C0, dnVDCC_C1, dnVDCC_C2, dnVDCC_C3, dnVDCC_C4, dnVDCC_O]
            self.dX[0] += dCa

        ### Calcium Buffers
        if 'calbindin' in self.models:
            dCa +=   + cbMoff*(cbH2M1 + cbH1M1 + cbH0M1 + 2*(cbH2M2 + cbH1M2 + cbH0M2)) \
                     + cbHoff*(cbH2M1 + cbH1M1 + cbH0M1 + 2*(cbH2M2 + cbH2M1 + cbH0M0)) \
                  - (+ cbMon*(cbH2M1 + cbH1M1 + cbH0M1  + 2*(cbH2M0 + cbH1M0 + cbH0M0)) \
                     + cbHon*(cbH2M1 + cbH1M1 + cbH0M1  + 2*(cbH0M2 + cbH0M1 + cbH0M0)) )*Ca

            dcbH0M0 = - 2*(cbMon+cbHon)*cbH0M0*Ca + cbMoff*cbH0M1 + cbHoff*cbH1M0

            dcbH0M1 = + (2*cbMon*cbH0M0 - (cbMon + 2*cbHon)*cbH0M1)*Ca \
                      + cbMoff*(2*cbH0M2 - cbH0M1) + cbHoff*cbH1M1

            dcbH0M2 = + cbMon*cbH0M1*Ca - 2*(cbMoff + cbHon*Ca)*cbH0M2 + cbHoff*cbH1M2

            dcbH1M0 = + (2*cbHon*cbH0M0 - (cbHon + 2*cbMon)*cbH1M0)*Ca\
                      + cbMoff*cbH1M1   + cbHoff*(2*cbH2M0 - cbH1M0)

            dcbH1M1 = - ((cbMon + cbHon)*Ca + (cbHoff + cbMoff))*cbH1M1\
                      + 2*((cbMoff*cbH1M2 + cbHoff*cbH2M1) + (cbMon*cbH1M0 + cbHon*cbH0M1)*Ca)

            dcbH1M2 = + (cbMon*cbH1M1 + 2*cbHon*cbH0M2 - cbHon*cbH1M2)*Ca\
                      - (2*cbMoff + cbHoff)*cbH1M2 + 2*cbHoff*cbH2M2

            dcbH2M0 = + cbMoff*cbH2M1 + cbHon*cbH1M0*Ca - 2*(cbHoff + cbMon*Ca)*cbH2M0

            dcbH2M1 = + (cbHon*cbH1M1 - cbMon*cbH2M1 + 2*cbMon*cbH2M0)*Ca\
                      + 2*(cbMoff*cbH2M2 - cbHoff*cbH2M1) - cbMoff*cbH2M1

            dcbH2M2 = + (cbMon*cbH2M1 + cbHon*cbH1M2)*Ca - 2*(cbHoff + cbMoff)*cbH2M2

            self.dX += [dcbH0M0, dcbH0M1, dcbH0M2, dcbH1M0, dcbH1M1, dcbH1M2, dcbH2M0, dcbH2M1, dcbH2M2]
            self.dX[0] += dCa

        ### Calcium Sensors
        if 'AZ' in self.models:
            #'''
            dCa += + (ab - af*Ca)*(AZ01 + AZ11 + AZ21 + AZ31 + AZ41 + AZ51)\
                   + 2*ab*b*(AZ02 + AZ12 + AZ22 + AZ32 + AZ42 + AZ52)\
                   - 2*af*(AZ00 + AZ10 + AZ20 + AZ30 + AZ40 + AZ50)*Ca\
                   + sb*((AZ10 + AZ11 + AZ12) + 2*b*(AZ20 + AZ21 + AZ22)\
                     + 3*b**2*(AZ30 + AZ31 + AZ32) + 4*b**3*(AZ40 + AZ41 + AZ42)\
                     + 5*b**4*(AZ50 + AZ51 + AZ52))\
                   - sf*((AZ40 + AZ42 + AZ41) + 2*(AZ30 + AZ31 + AZ32)\
                     + 3*(AZ20 + AZ21 + AZ22) + 4*(AZ10 + AZ12 + AZ11)\
                     + 5*(AZ00 + AZ01 + AZ02))*Ca

            dAZ00 =  + sb*AZ10 + ab*AZ01 - (5*sf + 2*af)*AZ00*Ca
            dAZ10 =  + 2*sb*b*AZ20 + 5*sf*AZ00*Ca + ab*AZ11\
                      - ((2*af + 4*sf)*Ca + sb)*AZ10
            dAZ20 =  + 3*sb*b**2*AZ30 + 4*sf*AZ10*Ca + ab*AZ21\
                      - ((2*af + 3*sf)*Ca + 2*sb*b)*AZ20
            dAZ30 =  + 4*sb*b**3*AZ40 + 3*sf*AZ20*Ca + ab*AZ31\
                      - (2*(af + sf)*Ca + 3*sb*b**2)*AZ30
            dAZ40 =  + 5*sb*b**4*AZ50 + 2*sf*AZ30*Ca + ab*AZ41\
                      - ((2*af + sf)*Ca + 4*sb*b**3)*AZ40
            dAZ50 =  + ab*AZ51 + (sf*AZ40 - 2*af*AZ50)*Ca - 5*sb*b**4*AZ50

            dAZ01 =  + 2*ab*b*AZ02 + sb*AZ11 + 2*af*AZ00*Ca\
                      - ((af + 5*sf)*Ca + ab)*AZ01
            dAZ11 =  + 2*(ab*AZ12 + sb*AZ21)*b + (5*sf*AZ01 + 2*af*AZ10)*Ca\
                      - ((af + 4*sf)*Ca + ab + sb)*AZ11
            dAZ21 =  + 2*ab*b*AZ22 + 3*sb*b**2*AZ31 + (4*sf*AZ11 + 2*af*AZ20)*Ca\
                      - ((af + 3*sf)*Ca + ab*AZ21 + 2*sb*b)*AZ21
            dAZ31 =  + 2*ab*b*AZ32 + 4*sb*b**3*AZ41 + (2*af*AZ30 + 3*sf*AZ21)*Ca\
                      - ((af + + 2*sf)*Ca + ab + 2*sf*Ca + 3*sb*b**2)*AZ31
            dAZ41 =  + 2*ab*b*AZ42 + 5*sb*b**4*AZ51 + 2*(sf*AZ31 + af*AZ40)*Ca\
                      - ((af + sf)*Ca + ab + 4*sb*b**3)*AZ41
            dAZ51 =  + 2*ab*b*AZ52 + (2*af*AZ50 + sf*AZ41)*Ca\
                      - (af*Ca + ab + 5*sb*b**4)*AZ51

            dAZ02 =  + af*AZ01*Ca + sb*AZ12 - (2*ab*b + 5*sf*Ca)*AZ02
            dAZ12 =  + (af*AZ11 + 5*sf*AZ02)*Ca + 2*sb*b*AZ22\
                      - (2*ab*b + 4*sf*Ca + sb)*AZ12
            dAZ22 =  + (af*AZ21 + 4*sf*AZ12)*Ca + 3*sb*b**2*AZ32\
                      - (2*(ab + sb)*b + 3*sf*Ca)*AZ22
            dAZ32 =  + 4*sb*b**3*AZ42 + (af*AZ31 + 3*sf*AZ22)*Ca\
                      - (2*(ab*b + sf*Ca) + 3*sb*b**2)*AZ32
            dAZ42 =  + (af*AZ41 + 2*sf*AZ32)*Ca + 5*sb*b**4*AZ52\
                      - (2*ab*b + sf*Ca + 4*sb*b**3)*AZ42
            dAZ52 =  + af*AZ51*Ca + sf*AZ42*Ca - 2*ab*b*AZ52 - 5*sb*b**4*AZ52


            self.dX += [dAZ00, dAZ10, dAZ20, dAZ30, dAZ40, dAZ50, dAZ01, dAZ11, dAZ21,\
                        dAZ31, dAZ41, dAZ51, dAZ02, dAZ12, dAZ22, dAZ32, dAZ42, dAZ52]
            self.dX[0] += dCa

        return  self.dX
