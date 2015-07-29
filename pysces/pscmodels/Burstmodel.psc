# To simulate install PySCeS 0.8.0 or newer (http://pysces.sf.net) 
# and the latest version of StochPy (http://stompy.sf.net)
 
# Keywords
Description: Burstmodel as used in Dobrzynski, M. and Bruggeman, F. J. (2009) Proc Natl Acad Sci U S A, 106(8), 2583 - 2588.
ModelType: Stochastic
Modelname: Burstmodel
Output_In_Conc: False
Species_In_Conc: True
 
# GlobalUnitDefinitions
UnitVolume: litre, 1.0, 0, 1
UnitLength: metre, 1.0, 0, 1
UnitSubstance: mole, 1.0, 0, 1
UnitTime: second, 1.0, 0, 1
UnitArea: metre, 1.0, 0, 2
 
# Reactions
R1:
    ONstate > OFFstate
    koff*ONstate

R2:
    OFFstate > ONstate
    kon*OFFstate

R3:
    $pool > mRNA
    ksyn*ONstate
# R3 has modifier(s): ONstate  

R4:
    mRNA > $pool
    kdeg*mRNA
 
# Fixed species
 
# InitPar
kon  = 0.05
koff = 0.80
kdeg = 2.5
ksyn = 80

# InitVar
ONstate = 0
OFFstate = 1
mRNA = 0

 

