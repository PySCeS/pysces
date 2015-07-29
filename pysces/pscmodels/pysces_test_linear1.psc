# PySCeS test input file
# Simple linear pathway (2004)

FIX: x0 x3

R1:
    x0 = s0
    k1*x0 - k2*s0

R2:
    s0 = s1
    k3*s0 - k4*s1

R3:
    s1 = s2
    k5*s1 - k6*s2

R4:
    s2 = x3
    k7*s2 - k8*x3

# InitExt
x0 = 10.0
x3 = 1.0

# InitPar
k1 = 10.0
k2 = 1.0
k3 = 5.0
k4 = 1.0
k5 = 3.0
k6 = 1.0
k7 = 2.0
k8 = 1.0

# InitVar
s0 = 1.0
s1 = 1.0
s2 = 1.0
