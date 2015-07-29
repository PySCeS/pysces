# PySCeS test input file
# Branched pathway (2004)

FIX: x0 x5 x6

R0:
    x0 = s1
    Vf0*(x0 - s1/Keq0)/(x0 + KS0*(1 + s1/KP0))


R1:
    s1 = s2
    Vf1*(s1 - s2/Keq1)/(s1 + KS1*(1 + s2/KP1))

R2:
    s2 = s3
    Vf2*(s2 - s3/Keq2)/(s2 + KS2*(1 + s3/KP2))

R3:
    s2 = s4
    Vf3*(s2 - s4/Keq3)/(s2 + KS3*(1 + s4/KP3))


R4:
    s3 = x5
    Vf4*(s3 - x5/Keq4)/(s3 + KS4*(1 + x5/KP4))

R5:
    s4 = x6
    Vf5*(s4 - x6/Keq5)/(s4 + KS5*(1 + x6/KP5))

# InitExt
x0 = 10.0
x5 = 1.0
x6 = 1.0

# InitPar
Vf0 = 10.0
Vf1 = 10.0
Vf2 = 10.0
Vf3 = 10.0
Vf4 = 10.0
Vf5 = 10.0

Keq0 = 10.0
Keq1 = 10.0
Keq2 = 10.0
Keq3 = 10.0
Keq4 = 10.0
Keq5 = 10.0

KS0 = 5.0
KS1 = 5.0
KS2 = 5.0
KS3 = 5.0
KS4 = 5.0
KS5 = 5.0

KP0 = 1.0
KP1 = 1.0
KP2 = 1.0
KP3 = 1.0
KP4 = 1.0
KP5 = 1.0

# InitVar
s1 = 1.0
s2 = 1.0
s3 = 1.0
s4 = 1.0
