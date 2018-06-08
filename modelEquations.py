from parameters import *

#print cmpts
class mfb:
    def __init__(self, models, name, dim):
        self.name = name
        self.dim  = dim
        self.vol  = reduce(lambda x, y: x*y, dim[3:])
        #print nbrs, dim, self.vol

        # Indexing of the compartment variables
        i = 0
        self.idx = od()
        self.X0 = []
        self.models = models
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

        if 'VDCC' in models:
            if models['VDCC'] == []: self.X0 += initVal['VDCC']
            else: self.X0 += models['VDCC']
            self.idx.update(od(zip(['VDCC_C0', 'VDCC_C1', 'VDCC_C2', 'VDCC_C3', 'VDCC_O'],
                                    range(i,i+5))))
            i += len(initVal['VDCC'])
            self.V = 0

        if 'calbindin' in models:
            if models['calbindin'] == []: self.X0 += initVal['calbindin']
            else: self.X0 += models['calbindin']
            self.idx.update(
                        od([('cbH0M0', i)  , ('cbH0M1', i+1), ('cbH0M2', i+2),
                            ('cbH1M0', i+3), ('cbH1M1', i+4), ('cbH1M2', i+5),
                            ('cbH2M0', i+6), ('cbH2M1', i+7), ('cbH2M2', i+8)]))
            i += len(initVal['calbindin'])

        if 'caSensor' in models:
            if models['caSensor'] == []: self.X0 += initVal['caSensor']
            else: self.X0 += models['caSensor']
            self.idx.update(
                    od([('CaS00', i+0 ), ('CaS10', i+1 ), ('CaS20', i+2 ),
                        ('CaS30', i+3 ), ('CaS40', i+4 ), ('CaS50', i+5 ),
                        ('CaS01', i+6 ), ('CaS11', i+7 ), ('CaS21', i+8 ),
                        ('CaS31', i+9 ), ('CaS41', i+10), ('CaS51', i+11),
                        ('CaS02', i+12), ('CaS12', i+13), ('CaS22', i+14),
                        ('CaS32', i+15), ('CaS42', i+16), ('CaS52', i+17)]))
            i += len(initVal['caSensor'])

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

        if 'VDCC' in self.models:
            j = len(initVal['VDCC'])
            VDCC_C0, VDCC_C1, VDCC_C2, VDCC_C3, VDCC_O = X[i:i+j]
            i += j

        if 'caSensor' in self.models:
            j = len(initVal['caSensor'])
            CaS00, CaS10, CaS20, CaS30, CaS40, CaS50, CaS01, CaS11, CaS21, \
            CaS31, CaS41, CaS51, CaS02, CaS12, CaS22, CaS32, CaS42, CaS52 = X[i:i+j]
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

        ### VDCC
        if 'VDCC' in self.models:
            dCa += 19.3*self.V*(0.3993 - np.exp(-self.V/80.36))/(1 - np.exp(self.V/80.36))*VDCC_O
            dVDCC_C0 = + b1(self.V)*VDCC_C1 - a1(self.V)*VDCC_C0
            dVDCC_C1 = + a1(self.V)*VDCC_C0 + b2(self.V)*VDCC_C2 \
                       - (b1(self.V) + a2(self.V))*VDCC_C1
            dVDCC_C2 = + a2(self.V)*VDCC_C1 + b3(self.V)*VDCC_C3 \
                       - (b2(self.V) + a3(self.V))*VDCC_C2
            dVDCC_C3 = + a3(self.V)*VDCC_C2 + b4(self.V)*VDCC_O \
                       - (b3(self.V) + a4(self.V))*VDCC_C3
            dVDCC_O  = + a4(self.V)*VDCC_C3 - b4(self.V)*VDCC_O

            self.dX += [dVDCC_C0, dVDCC_C1, dVDCC_C2, dVDCC_C3, dVDCC_O]
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
        if 'caSensor' in self.models:
            #'''
            dCa += + (ab - af*Ca)*(CaS01 + CaS11 + CaS21 + CaS31 + CaS41 + CaS51)\
                   + 2*ab*b*(CaS02 + CaS12 + CaS22 + CaS32 + CaS42 + CaS52)\
                   - 2*af*(CaS00 + CaS10 + CaS20 + CaS30 + CaS40 + CaS50)*Ca\
                   + sb*((CaS10 + CaS11 + CaS12) + 2*b*(CaS20 + CaS21 + CaS22)\
                     + 3*b**2*(CaS30 + CaS31 + CaS32) + 4*b**3*(CaS40 + CaS41 + CaS42)\
                     + 5*b**4*(CaS50 + CaS51 + CaS52))\
                   - sf*((CaS40 + CaS42 + CaS41) + 2*(CaS30 + CaS31 + CaS32)\
                     + 3*(CaS20 + CaS21 + CaS22) + 4*(CaS10 + CaS12 + CaS11)\
                     + 5*(CaS00 + CaS01 + CaS02))*Ca

            dCaS00 =  + sb*CaS10 + ab*CaS01 - (5*sf + 2*af)*CaS00*Ca
            dCaS10 =  + 2*sb*b*CaS20 + 5*sf*CaS00*Ca + ab*CaS11\
                      - ((2*af + 4*sf)*Ca + sb)*CaS10
            dCaS20 =  + 3*sb*b**2*CaS30 + 4*sf*CaS10*Ca + ab*CaS21\
                      - ((2*af + 3*sf)*Ca + 2*sb*b)*CaS20
            dCaS30 =  + 4*sb*b**3*CaS40 + 3*sf*CaS20*Ca + ab*CaS31\
                      - (2*(af + sf)*Ca + 3*sb*b**2)*CaS30
            dCaS40 =  + 5*sb*b**4*CaS50 + 2*sf*CaS30*Ca + ab*CaS41\
                      - ((2*af + sf)*Ca + 4*sb*b**3)*CaS40
            dCaS50 =  + ab*CaS51 + (sf*CaS40 - 2*af*CaS50)*Ca - 5*sb*b**4*CaS50

            dCaS01 =  + 2*ab*b*CaS02 + sb*CaS11 + 2*af*CaS00*Ca\
                      - ((af + 5*sf)*Ca + ab)*CaS01
            dCaS11 =  + 2*(ab*CaS12 + sb*CaS21)*b + (5*sf*CaS01 + 2*af*CaS10)*Ca\
                      - ((af + 4*sf)*Ca + ab + sb)*CaS11
            dCaS21 =  + 2*ab*b*CaS22 + 3*sb*b**2*CaS31 + (4*sf*CaS11 + 2*af*CaS20)*Ca\
                      - ((af + 3*sf)*Ca + ab*CaS21 + 2*sb*b)*CaS21
            dCaS31 =  + 2*ab*b*CaS32 + 4*sb*b**3*CaS41 + (2*af*CaS30 + 3*sf*CaS21)*Ca\
                      - ((af + + 2*sf)*Ca + ab + 2*sf*Ca + 3*sb*b**2)*CaS31
            dCaS41 =  + 2*ab*b*CaS42 + 5*sb*b**4*CaS51 + 2*(sf*CaS31 + af*CaS40)*Ca\
                      - ((af + sf)*Ca + ab + 4*sb*b**3)*CaS41
            dCaS51 =  + 2*ab*b*CaS52 + (2*af*CaS50 + sf*CaS41)*Ca\
                      - (af*Ca + ab + 5*sb*b**4)*CaS51

            dCaS02 =  + af*CaS01*Ca + sb*CaS12 - (2*ab*b + 5*sf*Ca)*CaS02
            dCaS12 =  + (af*CaS11 + 5*sf*CaS02)*Ca + 2*sb*b*CaS22\
                      - (2*ab*b + 4*sf*Ca + sb)*CaS12
            dCaS22 =  + (af*CaS21 + 4*sf*CaS12)*Ca + 3*sb*b**2*CaS32\
                      - (2*(ab + sb)*b + 3*sf*Ca)*CaS22
            dCaS32 =  + 4*sb*b**3*CaS42 + (af*CaS31 + 3*sf*CaS22)*Ca\
                      - (2*(ab*b + sf*Ca) + 3*sb*b**2)*CaS32
            dCaS42 =  + (af*CaS41 + 2*sf*CaS32)*Ca + 5*sb*b**4*CaS52\
                      - (2*ab*b + sf*Ca + 4*sb*b**3)*CaS42
            dCaS52 =  + af*CaS51*Ca + sf*CaS42*Ca - 2*ab*b*CaS52 - 5*sb*b**4*CaS52


            self.dX += [dCaS00, dCaS10, dCaS20, dCaS30, dCaS40, dCaS50, dCaS01, dCaS11, dCaS21,\
                        dCaS31, dCaS41, dCaS51, dCaS02, dCaS12, dCaS22, dCaS32, dCaS42, dCaS52]
            self.dX[0] += dCa

        return  self.dX
