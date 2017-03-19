# Tutorial: Kinetic Model
# lin4

FIX: X0 X4 S3

R1:
    X0 = S1
    (Vf1*(X0/Ks1) - Vr1*(S1/Kp1))/(1 + X0/Ks1 + S1/Kp1)

R2:
    S1 = S2
    (Vf2*(S1/Ks2) - Vr2*(S2/Kp2))/(1 + S1/Ks2 + S2/Kp2)

R3:
    S2 = S3
    random.gauss(S2, 0.1)

R4:
    S3 = X4
    S3*lognormal

!F lognormal = random.lognormvariate(2,1)*0.1

#InitEx
X0 = 10.0
X4 = 1.0

#InitPar
lognormal = 0

Vf1 = 10.0
Vr1 = 1.0
Ks1 = 1.0
Kp1 = 1.0 

Vf2 = 10.0
Vr2 = 1.0
Ks2 = 1.0 
Kp2 = 1.0 

Vf3 = 1.0
Ks3 = 1.0

#InitVar
S1 = 1
S2 = 1
S3 = 1
