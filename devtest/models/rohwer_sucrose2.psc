# Title Sucrose metabolism in sugar cane
# New version for resubmission of paper
# based on "sucrose2.cmd"
# Changed rate equations of R3 and R4 to have
# the same denominator, as it is the same enzyme

# Original model as published in Rohwer & Botha (2001) Biochem J 358:437-445

FIX: Fru_ex Glc_ex ATP ADP UDP phos glycolysis Suc_vac


R1:
     Fru_ex > Fru                # fructose uptake
     Vmax1*Fru_ex/(Km1_Fru_ex*(1 + Fru/Ki1_Fru) + Fru_ex)
                     # irr M-M with competitive inh by Fru

R2:
     Glc_ex > Glc                # glucose uptake
     Vmax2*Glc_ex/(Km2_Glc_ex*(1 + Glc/Ki2_Glc) + Glc_ex)
                     # irr M-M with competitive inh by Glc

R3:
     Glc + ATP = Hex_P + ADP    # glucokinase
# new rate equation below
     Vmax3*(Glc/Km3_Glc)*(ATP/Km3_ATP)/((1 + ATP/Km3_ATP) *
     (1 + Glc/Km3_Glc + Fru/Km4_Fru + G6P/Ki3_G6P + F6P/Ki4_F6P))
#     Vmax3*(Glc/Km3_Glc)*(ATP/Km3_ATP)/(1 + Glc/Km3_Glc + ATP/Km3_ATP
#     + Glc*ATP/(Km3_Glc*Km3_ATP) + G6P/Ki3_G6P + Fru/Km4_Fru + F6P/Ki4_F6P)
                     # random bireacant with competitive inh by
                     # G6P (product) and Fru

R4:
    Fru + ATP = Hex_P + ADP    # fructokinase 1 (Fru and Glc phosphorylating)
# new rate equation below
     Vmax4*(Fru/Km4_Fru)*(ATP/Km4_ATP)/((1 + ATP/Km4_ATP ) *
     (1 + Glc/Km3_Glc + Fru/Km4_Fru + G6P/Ki3_G6P + F6P/Ki4_F6P))
#     Vmax4*(Fru/Km4_Fru)*(ATP/Km4_ATP)/(1 + Fru/Km4_Fru + ATP/Km4_ATP
#     + Fru*ATP/(Km4_Fru*Km4_ATP) + F6P/Ki4_F6P + Glc/Km3_Glc + G6P/Ki3_G6P)
                     # random bireacant with competitive inh by
                     # F6P (product) and Glc

R5:
     Fru + ATP = Hex_P + ADP    # fructokinase 2 (only Fru phosphorylating)
#     Vmax5/(1 + F6P/Ki5_F6P + Fru/Ki5_Fru)*(Fru/Km5_Fru)*(ATP/Km5_ATP)/(1 +
     Vmax5/(1 + Fru/Ki5_Fru)*(Fru/Km5_Fru)*(ATP/Km5_ATP)/(1 +
     Fru/Km5_Fru + ATP/Km5_ATP + Fru*ATP/(Km5_Fru*Km5_ATP) + ADP/Ki5_ADP)
                     # random bireacant competitive inh by ADP
                     # non-competitive inh by F6P and Fru (delete F6P inh here!)

R6:
     {2}Hex_P = Suc6P + UDP     # SPS
     Vmax6f*(F6P*UDPGlc - Suc6P*UDP/Keq6)/(F6P*UDPGlc*(1+Suc6P/Ki6_Suc6P) +
     Km6_F6P*(1+Pi/Ki6_Pi)*(UDPGlc+Ki6_UDPGlc) + Km6_UDPGlc*F6P +
     Vmax6f/(Vmax6r*Keq6)*(Km6_UDP*Suc6P*(1+UDPGlc/Ki6_UDPGlc) +
     UDP*(Km6_Suc6P*(1+Km6_UDPGlc*F6P/(Ki6_UDPGlc*Km6_F6P*(1+Pi/Ki6_Pi))) +
     Suc6P*(1+F6P/Ki6_F6P))))
                     # ordered bi-bi
                     # Pi comp inh with respect to F6P

R7:
     Suc6P > Suc + phos              # SP
     Vmax7*Suc6P/(Km7_Suc6P + Suc6P)
                     # irr M-M

R8:
     Hex_P + Fru = Suc + UDP         # SS
     -Vmax8f*(Suc*UDP - Fru*UDPGlc/Keq8)/(Suc*UDP*(1+Fru/Ki8_Fru) +
     Km8_Suc*(UDP+Ki8_UDP) + Km8_UDP*Suc +
     Vmax8f/(Vmax8r*Keq8)*(Km8_UDPGlc*Fru*(1+UDP/Ki8_UDP) +
     UDPGlc*(Km8_Fru*(1+Km8_UDP*Suc/(Ki8_UDP*Km8_Suc)) +
     Fru*(1+Suc/Ki8_Suc))))
                     # ordered bi-bi

R9:
     Suc > Glc + Fru                 # INV
     Vmax9/(1+Glc/Ki9_Glc)*Suc/(Km9_Suc*(1+Fru/Ki9_Fru) + Suc)
                     # irr M-M with comp inh by Fru
                     # and non-comp by Glc

R10:
    Hex_P > glycolysis             # glycolytic drain
    Vmax10*F6P/(Km10_F6P + F6P)
                     # irr M-M

R11:
    Suc > Suc_vac                 # Suc import into vacuole (active)
    Vmax11*Suc/(Km11_Suc + Suc)
                     # irr M-M

#------------------kinetic constants---------------------------

# NB all rates in mM/min, assuming 1 g FW = 700 ul internal volume
# 90% of this is vacuole, so that 1 g FW = 70 ul cytosolic vol

# all concentrations in mM 

# R1 
Vmax1         = 0.286 
Km1_Fru_ex    = 0.2 
Ki1_Fru       = 1.0     # estimate, will allow high accumulation (e.g.
                       # 30 mM in internode 5) 

# R2 
Vmax2         = 0.286  # same as Fru 
Km2_Glc_ex    = 0.2 
Ki2_Glc       = 1.0 

# R3 
Vmax3         = 0.197 
Km3_Glc       = 0.07   # range 35 - 100 uM 
Km3_ATP       = 0.25   # range 90 - 500 uM 
Ki3_G6P       = 0.1 

# R4 
Vmax4         = 0.197  # same Vmax for Glc and Fru 
Km4_Fru       = 10.0 
Km4_ATP       = 0.25   # same enzyme 
Ki4_F6P       = 10.0     # estimate, same as Km_Fru 

# R5 
Vmax5         = 0.164 
Km5_Fru       = 0.1 
Km5_ATP       = 0.085 
Ki5_ADP       = 2.0 
# Ki5_F6P       = 1.3 
Ki5_Fru       = 12.0     # range 5 - 21 mM 

# R6 
Keq6          = 10.0    # ##double check!!delG = -5.7 kJ/mol 
Vmax6f        = 0.379 
Vmax6r        = 0.2    # estimate 
Km6_F6P       = 0.6 
Km6_UDPGlc    = 1.8 
Km6_UDP       = 0.3    # estimate, measured [] = Km
Km6_Suc6P     = 0.1    # estimate, same order of magnitued as SP 
Ki6_Pi        = 3.0      # range 1 - 5 mM 
Ki6_F6P       = 0.4    # estimate, slightly less than Km 
Ki6_UDPGlc    = 1.4    # estimate, slightly less than Km 
Ki6_Suc6P     = 0.07   # estimate, slightly less than Km 

# R7 
Vmax7         = 0.5    # estimate, higher than SPS 
Km7_Suc6P     = 0.1    # range 45 - 150 uM 

# R8 
Keq8          = 5.0      # delG = -4 kJ/mol 
Vmax8f        = 0.677 
Vmax8r        = 0.3    # estimate, check later 
Km8_Suc       = 50     # one report 87, other 5-50 
Km8_UDP       = 0.3    # range 212-390 uM 
Km8_UDPGlc    = 0.3    # range 0.06-1.7?? (Avigad 1982) 
Km8_Fru       = 4.0      # range 1.4-6.2 (Avigad 1982) 
Ki8_Suc       = 40.0     # estimate 
Ki8_Fru       = 4.0      # range 2.48-17.8 
Ki8_UDP       = 0.3    # same as Km 

# R9 
Vmax9         = 0.372 
Km9_Suc       = 10.0 
Ki9_Fru       = 15.0     # range 8.6-30 
Ki9_Glc       = 15.0 

# R10 
Vmax10        = 0.1    # tuning handle 
Km10_F6P      = 0.2    # Km of PFP 

# R11 
Vmax11        = 1.0      # tuning handle 
Km11_Suc      = 100.0    # thumbsuck 

#--------------terminal fixed metabolites---------------------
# internode 5 
Fru_ex        = 5.0 
Glc_ex        = 5.0 
ATP           = 1.0 
ADP           = 0.2 
UDP           = 0.2   # assume ATP=UTP and ADP=UDP 
Pi            = 5.1 
phos          = 5.1 
glycolysis    = 0.0 
Suc_vac       = 0.0 

#--------------variable metabolites---------------------------
Fru           = 1.0 
Glc           = 1.0 
Hex_P         = 1.0 
UDPGlc        = 1.0 
G6P           = 1.0 
F6P           = 1.0 
Suc6P         = 1.0 
Suc           = 1.0 



#-------------- FORCING FUNCTIONS FOR HEXOSE PHOSPHATE BLOCK ----------------
# specified with !F in input file
!F self.UDPGlc = 0.8231*self.Hex_P 
!F self.G1P    = 0.0064*self.Hex_P 
!F self.G6P    = 0.1130*self.Hex_P 
!F self.F6P    = 0.0575*self.Hex_P 

# for Vmax8f scan
!F self.Vmax8r = self.Vmax8f/2.2567

# for Vmax6f scan
!F self.Vmax6r = self.Vmax6f/1.895

#-------------- INITIALISATION FUNCTIONS TO FORCE NUMERIC DERIVATION OF ELASTICITIES ----------------
# specified with !I in input file
#!I self.mode_elas_deriv=1
