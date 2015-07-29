FIX: S DEM

R1:
    S = A
    V1*(S/S05)*(1 - A/(S*Keq1))*((S/S05)+(A/A05))**(h-1)/((S/S05 + A/A05)**h + (1 + (X/X05)**h)/(1 + alpha*(X/X05)**h))

R2:
    A = B
    V2*(A - B/Keq2)/(A + K2A*(1 + B/K2B))

R3:
    B = X
    V3*(B - X/Keq3)/(B + K3B*(1 + X/K3X))

R4:
    X = DEM
    V4*X/(X + K4X)

Keq1 = 400.0
V1 = 200.0
S05 = 1.0
A05 = 10000 #p05
X05 = 1
h = 4.0
alpha = 0.01

V2 = 100.00
Keq2 = 10.0
K2A = 1.0
K2B = 1.0

V3 = 100.0
Keq3 = 10.0
K3B = 1.0
K3X = 1.0

V4 = 100.0
K4X = 0.1


S = 1.0
A = 1.0
B = 1.0
X = 1.0
DEM = 1.0
