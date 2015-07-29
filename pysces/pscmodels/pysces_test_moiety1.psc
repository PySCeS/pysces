# PySCeS test input file
# Moiety conserved pathway (2004)

FIX: x0 x1 x2 x3 x4 x5

R1:
x2 + s4 = s9 + s5
k1*x2*s4 - k2*s9*s5

R2:
s9 + s3 = s7 + s2
k3*s9*s3 - k4*s7*s2

R3:
s7 + s0 = x3 + s1
k5*s7*s0 - k6*x3*s1

R4:
x5 + s2 = s6 + s3
k7*x5*s2 - k8*s6*s3

R5:
s6 + s5 = x4 + s4
k9*s6*s5 - k10*x4*s4

R6:
x1 + s1 = s8 + s0
k11*x1*s1 - k12*s8*s0

R7:
s8 = x0
k13*s8 - k14*x0

# InitPar
k1 = 10.0
k2 = 1.0
k3 = 5.0
k4 = 1.0
k5 = 15.0
k6 = 1.0
k7 = 10.0
k8 = 1.0
k9 = 15.0
k10 = 1.0
k11 = 10.0
k12 = 1.0
k13 = 15.0
k14 = 1.0

# InitExt
x0 = 5.0
x1 = 10.0
x2 = 10.0
x3 = 5.0
x4 = 5.0
x5 = 10.0

# InitVar
s0 = 1.0
s1 = 2.0
s2 = 3.0
s3 = 4.0
s4 = 5.0
s5 = 6.0
s6 = 7.0
s7 = 8.0
s8 = 9.0
s9 = 10.0
