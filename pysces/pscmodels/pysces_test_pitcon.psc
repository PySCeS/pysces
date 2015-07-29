# A three enzyme linear system with feedback
# which exhibits multiple steady-state solutions
# Brett Olivier, Johann Rohwer and Jannie Hofmeyr

FIX: S P

R1:
    S = A
    V1*(S/S05)*(1 - A/(S*Keq1))*((S/S05)+(A/A05))**(h-1)/
    ((S/S05 + A/A05)**h + (1 + (P/P05)**h)/(1 + alpha*(P/P05)**h))

R2:
    A = B
    V2*(A - B/Keq2)/(A + K2A*(1 + B/K2B))

R3:
    B = P
    V3*(B - P/Keq3)/(B + K3B*(1 + P/K3P))

Keq1 = 400.0
V1 = 200.0
S05 = 1.0
A05 = 10000.0
P05 = 1.0
h = 4.0
alpha = 0.01

V2 = 10000.0
Keq2 = 10.0
K2A = 1.0
K2B = 1.0

V3 = 10000.0
Keq3 = 10.0
K3B = 1.0
K3P = 1.0

S = 1.0
A = 10.0
B = 10.0
P = 1.0
