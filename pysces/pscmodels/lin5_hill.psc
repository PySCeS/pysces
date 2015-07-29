FIX: S P

R1:
    S = A
    V1*(S/S05)*(1 - A/(S*Keq1))*((S/S05)+(A/A05))**(h-1)/
    ((S/S05 + A/A05)**h + (1 + (C/C05)**h)/(1 + alpha*(C/C05)**h))

R2:
    A = B
    V2*(A - B/Keq2)/(A + K2A*(1 + B/K2B))

R3:
    B = C
    V3*(B - C/Keq3)/(B + K3B*(1 + C/K3C))
    
R4:
    C = D
    V4*(C - D/Keq4)/(C + K4C*(1 + D/K4D))
    
R5:
    D = P
    (V5*D)/(D + K5D)

Keq1 = 400.0
V1 = 200.0
S05 = 1.0
A05 = 10000.0
C05 = 1.0
h = 4.0
alpha = 0.01

V2 = 10000.0
Keq2 = 10.0
K2A = 1.0
K2B = 1.0

V3 = 10000.0
Keq3 = 10.0
K3B = 1.0
K3C = 1.0

V4 = 10.0
Keq4 = 10.0
K4C = 0.01
K4D = 1.0

V5 = 20.0
K5D = 1.0

S = 1.0
A = 1.0
B = 1.0
C = 1.0
D = 1.0
P = 1.0

#!U self.GoK1 = self.A/self.S/self.Keq1
#!U self.GoK2 = self.B/self.A/self.Keq2
#!U self.GoK3 = self.C/self.B/self.Keq3
#!U self.GoK4 = self.D/self.C/self.Keq4
#!U print 'GoK1 = ', self.GoK1
#!U print 'GoK2 = ', self.GoK2
#!U print 'GoK3 = ', self.GoK3
#!U print 'GoK4 = ', self.GoK4

