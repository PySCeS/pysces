.. _PySCeS-Inputfile:

The PySCeS Model Description Language
=====================================

 PySCeS: the Python Simulator for Cellular Systems is an
 extendable toolkit for the analysis and investigation of cellular
 systems. It is available for download from: http://pysces.sf.net

PySCeS uses an ASCII text based *input file* to describe a 
cellular system in terms of it's stoichiometry, kinetics, 
compartments and parameters. Input files may have any filename 
with the single restriction that, for cross platform 
compatibility, they must end with the extension *.psc*. In this 
document we describe the PySCeS Model Description Language  
(MDL) which has been updated and extended for the PySCeS 0.7.x 
release. 

 PySCeS is distributed under the PySCeS (BSD style) license and is made
 freely available as Open Source software. See LICENCE.txt for details.
 
We hope that you will enjoy using our software. If, however, 
you find any unexpected features (i.e. bugs) or have any 
suggestions on how we can improve PySCeS and specifically the 
PySCeS MDL please let us know. 

.. _PySCeS-Inputfile-Detailed:

Defining a PySCeS model
-----------------------

.. _PySCeS-Inputfile-Basic:

A kinetic model
~~~~~~~~~~~~~~~

The basic description of a kinetic model in the PySCeS MDL contains
the following information:

* whether any fixed (boundary) species are present
* the reaction network stoichiometry
* rate equations for each reaction step
* parameter and boundary species initial values
* the initial values of the variable species

Although it is in principle possible to define an ODE based model
without reactions or free species, for practical purposes PySCeS
requires a minimum of a single reaction. Once this information is
obtained it can be organised and written as a PySCeS input file.
While this list is the minimum information required for a PySCeS input
file the MDL allows the definition of advanced models that contain
compartments, global units, functions, rate and assignment rules.

.. _PySCeS-Inputfile-Detailed-Keywords:

Model keywords
~~~~~~~~~~~~~~

In PySCeS 0.7.x it is now possible to define keywords that
specify model information. Keywords have the general form

``<keyword>: <value>``

The *Modelname* (optional) keyword, containing only 
alphanumeric characters (or _), describes the model filename 
(typically used when the model is exported via the PySCeS 
interface module) while the *Description* keyword is a (short) 
single line model description. :: 

 Modelname: rohwer_sucrose1
 Description: Sucrose metabolism in sugar cane (Johann M. Rohwer)
 
Two keywords are available for use (optional) with models that have one or
more compartments defined. Both take a boolean (True/False) as
their value: 

 * *Species_In_Conc* specifies whether the species symbols used in the rate equations represent a concentration (True, default) or an amount (False).
 * *Output_In_Conc* tells PySCeS to output the results of numerical operations in concentrations (True, default) or in amounts (False).
 
::
 
  Species_In_Conc: True
  Output_In_Conc: False
 
More information on the effect these keywords have on the analysis of a model
can be found in the PySCeS Reference Manual. 
 
.. _PySCeS-Inputfile-Detailed-Units:

Global unit definition
~~~~~~~~~~~~~~~~~~~~~~

PySCeS 0.7 supports the (optional) definition of a set of 
global units. In doing so we have chosen to follow the general 
approach used in the Systems Biology Modelling Language (SBML 
L2V3) specification. The general definition of a PySCeS unit 
is: ```<UnitType>: <kind>, <multiplier>, <scale>, <exponent>``` 
where *kind* is a string describing the base unit (for SBML 
compatibility this should be an SI unit) e.g. mole, litre, 
second or metre. The base unit is modified by the multiplier, 
scale and index using the following relationship: 
*<multiplier> * (<kind> * 10**<scale>)**<index>*. The 
default unit definitions are:: 

 UnitSubstance: mole, 1, 0, 1
 UnitVolume: litre, 1, 0, 1
 UnitTime: second, 1, 0, 1
 UnitLength: metre, 1, 0, 1
 UnitArea: metre, 1, 0, 2

Please note that defining these values does not affect the 
numerical analysis of the model in any way. 

.. _PySCeS-Inputfile-Detailed-Names:

Symbol names and comments
~~~~~~~~~~~~~~~~~~~~~~~~~

Symbol names (i.e. reaction, species, compartment, function, 
rule and parameter names etc.) must start with either an 
underscore or letter and be followed by any combination of 
alpha-numeric characters or an underscore. Like all other 
elements of the input file names are case sensitive:: 

 R1
 _subA
 par1b
 ext_1

Explicit access to the "current" time in a time simulation is 
provided by the special symbol ``_TIME_``. This is useful in 
the definition of events and rules (see chapter on advanced 
model construction for more details). 

Comments can be placed anywhere in the input file in one of two 
ways, as single line comment starting with a *#* or as a 
multi-line triple quoted comment *"""<comment>"""*::

 # everything after this is ignored
 
 """
 This is a comment
 spread over a
 few lines.
 """



.. _PySCeS-Inputfile-Detailed-Compartments:

Compartment definition
~~~~~~~~~~~~~~~~~~~~~~

By default (as is the case in all PySCeS versions < 0.7) PySCeS 
assumes that the model exists in a single unit volume 
compartment. In this case it is **not** necessary to define a 
compartment and the ODE's therefore describe changes in 
concentration per time. However, if a compartment is defined, 
PySCeS assumes that the ODE's describe changes in substance amount per 
time. Doing this affects how the model is defined in the input 
file (especially with respect to the definitions of rate 
equations and species) and the user is **strongly** advised to 
read the Users Guide before building models in this way. The 
compartment definition is as follows ``Compartment: <name>, 
<size>, <dimensions>``, where *<name>* is the unique 
compartment id, *<size>* is the size of the compartment (i.e. 
length, volume or area) defined by the number of *<dimensions>* 
(e.g. 1,2,3):: 

 Compartment: Cell, 2.0, 3
 Compartment: Memb, 1.0, 2 

.. _PySCeS-Inputfile-Detailed-Functions:

Function definitions
~~~~~~~~~~~~~~~~~~~~

A new addition to the PySCeS MDL is the ability to define SBML 
styled functions. Simply put these are code substitutions that 
can be used in rate equation definitions to, for example, 
simplify the kinetic law. The general syntax for a function is 
``Function: <name>, <args> {<formula>}`` where *<name>* is the 
unique function id, *<arglist>* is one or more comma separated 
function arguments. The *<formula>* field, enclosed in curly 
brackets, may only make use of arguments listed in the 
*<arglist>* and therefore **cannot** reference model attributes 
directly. If this functionality is required a forcing function 
(assignment rule) may be what you are looking for. :: 

 Function: rmm_num, Vf, s, p, Keq {
 Vf*(s - p/Keq)
 }

 Function: rmm_den, s, p, Ks, Kp {
 s + Ks*(1.0 + p/Kp)
 }

The syntax for function definitions has been adapted from Frank 
Bergmann and Herbert Sauro's "Human Readable Model Definition 
Language" (Draft 1). 

.. _PySCeS-Inputfile-Detailed-Fixed:

Defining fixed species
~~~~~~~~~~~~~~~~~~~~~~

Boundary species, also known as fixed or external species, are 
a special class of parameter used when modelling biological 
systems. The PySCeS MDL fixed species are declared on a single 
line as ``FIX: <fixedlist>``. The *<fixedlist>* is a space 
separated list of symbol names which should be initialised like 
any other species or parameter::

 FIX: Fru_ex Glc_ex ATP ADP UDP phos glycolysis Suc_vac

If no fixed species are present in the model then this 
declaration should be omitted entirely. 

.. _PySCeS-Inputfile-Detailed-Reactions:

Reaction stoichiometry and rate equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reaction stoichiometry and rate equation are defined together
as a single reaction step. Each step in the system is defined as
having a name (identifier), a stoichiometry (substrates are
converted to products) and rate equation (the catalytic activity,
described in terms of species and parameters). All reaction
definitions should be separated by an empty line. The general format
of a reaction in a model with no compartments is::

 <name>: 
         <stoichiometry>
         <rate equation>

The *<name>* argument follows the syntax as discussed in a 
previous section, however, when more than one compartment has 
been defined it is important to locate the reaction in its 
specific compartment. This is done using the ``@`` operator:: 

 <name>@<compartment>: 
                       <stoichiometry>
                       <rate equation>

Where *<compartment>* is a valid compartment name. In either 
case this then followed either directly (or on the next line) 
by the reaction stoichiometry.  

Each *<stoichiometry>* argument is defined in terms of reaction 
substrates, appearing on the left hand side and products on the 
right hand side of an identifier which labels the reaction as 
either reversible (*=*) or irreversible (*>*). If required each 
reagent's stoichiometric coefficient (PySCeS accepts both 
integer and floating point) should be included in curly braces 
*{}* immediately preceding the reagent name. If these are 
omitted a coefficient of one is assumed:: 

 {2.0}Hex_P = Suc6P + UDP  # reversible reaction
 Fru_ex > Fru              # irreversible reaction
 species_5 > $pool         # a reaction to a sink

The PySCeS MDL also allows the use of the *$pool* token that 
represents a placeholder reagent for reactions that have no 
net substrate or product. Reversibility of a reaction is only 
used when exporting the model to other formats (such as SBML) 
and in the calculation of elementary modes. It does not affect 
the numerical evaluation of the rate equations in any way. 

Central to any reaction definition is the *<rate equation>* 
(SBML kinetic law). This should be written as valid Python 
expression and may fall across more than one line. Standard 
Python operators ``+ - * / **`` are supported (note the Python 
power e.g. *2^4* is written as *2\*\*4*). There is no shorthand 
for multiplication with a bracket so *-2(a+b)^h* would be written as 
*-2\*(a+b)\*\*h}* and normal operator precedence applies: 

 +--------+-------------------------+
 |  +, -  | addition, subtraction   |
 +--------+-------------------------+
 |  \*, / | multiplication, division|
 +--------+-------------------------+
 | +x,-x  | positive, negative      |
 +--------+-------------------------+
 |  \*\*  | exponentiation          |
 +--------+-------------------------+
 
Operator precedence increase from top to bottom and left to 
right (adapted from the Python Reference Manual). 

The PySCeS MDL parser has been developed to parse and translate different
styles of infix into Python/Numpy based expressions, the following
functions are supported in any mathematical expression:

 * log, log10, ln, abs
 * pow, exp, root, sqrt
 * sin, cos, tan, sinh, cosh, tanh
 * arccos, arccosh, arcsin, arcsinh, arctan, arctanh
 * floor, ceil, ceiling, piecewise
 * notanumber, pi, infinity, exponentiale

Logical operators are supported in rules, events etc but *not*
in rate equation definitions. The PySCeS parser understands
Python infix as well as libSBML and NumPy prefix notation.  

 * and or xor not
 * > gt(x,y) greater(x,y)
 * < lt(x,y) less(x,y)
 * >= ge(x,y) geq(x,y) greater_equal(x,y) 
 * <= le(x,y) leq(x,y) less_equal(x,y)
 * == eq(x,y) equal(x,y) 
 * != neq(x,y) not_equal(x,y)

Note that currently the MathML *delay and factorial* functions 
are not supported. Delay is handled by simply removing it from 
any expression, e.g. *delay(f(x), delay)* would be parsed as 
*f(x)*. Support for *piecewise* has been recently added 
to PySCeS and will be discussed in the *advanced features* section. 

A reaction definition when no compartments are defined::

 R5: Fru + ATP = Hex_P + ADP
     Fru/Ki5_Fru)*(Fru/Km5_Fru)*(ATP/Km5_ATP)/(1 +
     Vmax5/(1 + Fru/Ki5_Fru)*(Fru/Km5_Fru)*(ATP/Km5_ATP)/(1 +
     Fru/Km5_Fru + ATP/Km5_ATP + Fru*ATP/(Km5_Fru*Km5_ATP) +
     ADP/Ki5_ADP)

and using the previously defined functions::

 R6:
    A = B
    rmm_num(V2,A,B,Keq2)/rmm_den(A,B,K2A,K2B)

When compartments are defined note how now the reaction is now 
given a location and that because the ODE's formed from these 
reactions must be in changes in substance per time the rate 
equation is multiplied by its compartment size. In this 
particular example the species symbols represent concentrations 
(*Species_In_Conc: True*):: 

 R1@Cell:
     s1 = s2
     Cell*(Vf1*(s1 - s2/Keq1)/(s1 + KS1*(1 + s2/KP1)))

If *Species_In_Conc: True* the location of the species is 
defined when it is initialised and will be explained later in 
this manual. The following example shows the species symbols 
explicitly defined as amounts (*Species_In_Conc: False*):: 

 R4@Memb: s3 = s4
     Memb*(Vf4*((s3/Memb) - (s4/Cell)/Keq4)/((s3/Memb)
     + KS4*(1 + (s4/Cell)/KP4)))

Please note that at this time we are not certain if this form 
of rate equation is translatable into valid SBML in a way that is 
interoperable with other software. 

.. _PySCeS-Inputfile-Detailed-Initialisation:

Species and parameter initialisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The general form of any species (fixed, free) and parameter is 
simply:: 

 property = value

Initialisations can be written in any order anywhere in the 
input file but for human readability purposes these are usually 
placed after the reaction that uses them or grouped at the end 
of the input file. Both decimal and scientific notation is 
allowed with the following provisions that neither floating 
point *(1. )* nor scientific shorthand *(1.e-3)* syntax should 
be used, instead use the full form *(1.0e-3)*, *(0.001)* or 
*(1.0)*. 

Variable or free species are initialised differently depending 
on whether compartments are present in the model. While in 
essence the variables are set by the system parameters the 

Although the variable species concentrations are determined by 
the parameters of the system, their initial values are used in 
various places, calculating total moiety concentrations (if 
present), time simulation initial values (e.g. time=zero) and 
as initial guesses for the steady-state algorithms. If an empty 
initial species pool is required it is not recommended to 
initialise these values to zero (in order to prevent potential 
divide-by-zero errors) but rather to a small value (e.g. 
10**-8). 

For a model with no compartments these initial values assumed 
to be concentrations:: 

 NADH = 0.001
 ATP  = 2.3e-3
 sucrose = 1
 
In a model with compartments it is expected that the species 
are located in a compartment (even if *Species_In_Conc: False*) 
this is done useing the *@* symbol:: 

 s1@Memb = 0.01
 s2@Cell = 2.0e-4

A word of warning, the user is responsible for making sure that 
the units of the initialised species match those of the model. 
Please keep in mind that **all** species (and anything that 
depends on them) is defined in terms of the *Species_In_Conc* 
keyword. For example, if the preceding initialisations were for 
*R1* (see Reaction section) then they would be concentrations 
(as *Species_In_Conc: True*). However, in the next example, we 
are initialising species for *R4* and they are therefore in 
amounts (*Species_In_Conc: False*):: 

 s3@Memb = 1.0
 s4@Cell = 2.0

Fixed species are defined in a similar way and although 
technically a parameter, they should be given a location in 
compartmental models:: 

 # InitExt
 X0 = 10.0
 X4@Cell = 1.0

However, fixed species are true parameters in the sense that 
their associated compartment size does not affect their value 
when it changes size. If compartment size dependent behaviour 
is required an assignment or rate rule should be considered. 

Finally, the parameters should be initialised. PySCeS checks if 
a parameter is defined that is not present in the rate 
equations and if such parameter initialisations are detected a 
harmless warning is generated. If, on the other hand, an 
uninitialised parameter is detected a warning is generated and 
a value of 1.0 assigned:: 

 # InitPar
 Vf2 = 10.0
 Ks4 = 1.0

.. _PySCeS-Inputfile-Advanced:

Advanced model construction
---------------------------

.. _PySCeS-Inputfile-Advanced-Assignment:

Assignment rules
~~~~~~~~~~~~~~~~

Assignment rules or forcing functions are used to set the value 
of a model attribute before the ODE's are evaluated. This model 
attribute can either be a parameter used in the rate equations 
(this is traditionally used to describe an equilibrium block) a 
compartment or an arbitrary parameter (commonly used to define 
some sort of tracking function). Assignment rules can access 
other model attributes directly and have the generic form ``!F 
<par> = <formula>``. Where *<par>* is the parameter assigned 
the result of *<formula>*. Assignment rules can be defined 
anywhere in the input file:: 

 !F S_V_Ratio = Mem_Area/Vcyt
 !F sigma_test = sigma_P*Pmem + sigma_L*Lmem
 
These rules would set the value of *<par>* which whose value 
can be followed with using the simulation and steady state 
extra_data functionality. 

.. _PySCeS-Inputfile-Advanced-Raterule:

Rate rules
~~~~~~~~~~

PySCeS now includes support for rate rules which are 
essentially directly encoded ODE's which are evaluated after 
the ODE's defined by the model stoichiometry and rate 
equations. Unlike the SBML rate rule, PySCeS allows one to 
access a reaction symbol in the rate rules (this is 
automatically expanded when the model is exported to SBML). The 
general form of a rate rule is ``RateRule: <par> 
{<function>}``. Where *<name>* is the model attribute (e.g. 
compartment or parameter) whose rate of change is described by 
the *<formula>*. It may also be defined anywhere in the input 
file:: 

 RateRule: Mem_Area {
 (sigma_P)*(Mem_Area*k4*(P)) + (sigma_L)*(Mem_Area*k5*(L))
 }
 
 RateRule: Vcyt {(1.0/Co)*(R1()+(1-m1)*R2()+(1-m2)*R3()-R4()-R5())}

Remember to initialise any new parameters used in the rate rules.
 
.. _PySCeS-Inputfile-Advanced-Events:

Events
~~~~~~

Time dependant events may now be defined whose definition 
follows the event framework described in the SBML L2V1 
specification. The general form of an event is *Event: <name>, 
<trigger>, <delay> { <assignments> }*. As can be seen an event 
consists of essentially three parts, a conditional *<trigger>*, 
a set of one or more *<assignments>* and a *<delay>* between 
when the trigger is fired (and the assignments are evaluated) 
and the eventual assignment to the model. Assignments have the 
general form *<par> = <formula>*. Events have access to the 
"current" simulation time using the *_TIME_* symbol:: 

 Event: event1, _TIME_ > 10 and A > 150.0, 0 {
 V1 = V1*vfact
 V2 = V2*vfact
 }

The following event illustrates the use of a delay of ten time 
units as well as the prefix notation (used by libSBML) for the 
trigger (PySCeS understands both notations):: 

 Event: event2, geq(_TIME_, 15.0), 10 {
 V3 = V3*vfact2
 } 

*Note:* in order for PySCeS to handle events it is necessary to
have the PySundials installed

.. _PySCeS-Inputfile-Advanced-Piecewise:

Piecewise
~~~~~~~~~

Although technically an operator piecewise functions are 
sufficiently complicated to warrant their own section. A 
piecewise operator is essentially an *if, elif, ..., else* 
logical operator that can be used to conditionally "set" the 
value of some model attribute. Currently piecewise is supported 
in rule constructs and has not been tested directly in rate 
equation definitions. The piecewise function's most basic 
incarnation is `piecewise(<val1>, <cond>, <val2>)` which is evaluated as:: 

 if <cond>:
     return <val1>
 else:
     return <val2>

alternatively, `piecewise(<val1>, <cond1>, <val2>, <cond2>, <val3>, <cond3>)`::

 if <cond1>:
     return <val1>
 elif <cond2>:
     return <val1>
 elif <cond3>:
     return <val3>

or `piecewise(<val1>, <cond1>, <val2>, <cond2>, <val3>, <cond3>, <val4>)`::

 if <cond1>:
     return <val1>
 elif <cond2>:
     return <val2>
 elif <cond3>:
     return <val3>
 else:
     return <val4>

can also be used. A "real-life" example of an assignment rule 
with a piecewise function:: 

 !F Ca2plus=piecewise(0.1, lt(_TIME_,60), 0.1, gt(_TIME_,66.0115), 1)  

In principle there is no limit on the amount of conditional 
statements present in a piecewise function, the condition can 
be a compound statements *a or b and c* and may include the 
*_TIME_* symbol. 

Reagent placeholder
~~~~~~~~~~~~~~~~~~~

Some models contain reactions which are defined as only have substrates or
products::

 R1: A + B >
 
 R2: > C + D
 
The implication is that the relevant reagents appear or disappear
from or into a constant pool. Unfortunately the `PySCeS` parser does not accept
such an unbalanced reaction definition and requires these pools to be
represented as a ``$pool`` token::

 R1: A + B > $pool
 
 R2: $pool > C + D

``$pool`` is neither counted as a reagent nor does it ever appear in the
stoichiometry (think of it as dev/null) and no other $<str> tokens are allowed.


.. _PySCeS-Inputfile-Examples:

Example PySCeS input files
--------------------------

.. _PySCeS-Inputfile-Examples-Basic:

Basic model definition
~~~~~~~~~~~~~~~~~~~~~~

PySCeS test model *pysces_test_linear1.psc*:: 

 FIX: x0 x3

 R1: x0 = s0
     k1*x0 - k2*s0

 R2: s0 = s1
     k3*s0 - k4*s1

 R3: s1 = s2
     k5*s1 - k6*s2

 R4: s2 = x3
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

.. _PySCeS-Inputfile-Examples-Advanced:

Advanced example
~~~~~~~~~~~~~~~~

This model includes the use of *Compartments*, *KeyWords*, 
*Units* and *Rules*:: 

 Modelname: MWC_wholecell2c
 Description: Surovtsev whole cell model using J-HS Hofmeyr's framework

 Species_In_Conc: True
 Output_In_Conc: True

 # Global unit definition
 UnitVolume: litre, 1.0, -3, 1
 UnitSubstance: mole, 1.0, -6, 1
 UnitTime: second, 60, 0, 1

 # Compartment definition
 Compartment: Vcyt, 1.0, 3
 Compartment: Vout, 1.0, 3
 Compartment: Mem_Area, 5.15898, 2

 FIX: N 

 R1@Mem_Area: N = M
    Mem_Area*k1*(Pmem)*(N/Vout)

 R2@Vcyt: {244}M = P # m1
    Vcyt*k2*(M)

 R3@Vcyt: {42}M = L # m2
    Vcyt*k3*(M)*(P)**2

 R4@Mem_Area: P = Pmem
    Mem_Area*k4*(P)

 R5@Mem_Area: L = Lmem
    Mem_Area*k5*(L)

 # Rate rule definition
 RateRule: Vcyt {(1.0/Co)*(R1()+(1-m1)*R2()+(1-m2)*R3()-R4()-R5())}
 RateRule: Mem_Area {(sigma_P)*R4() + (sigma_L)*R5()}

 # Rate rule initialisation
 Co = 3.07e5 # uM p_env/(R*T)
 m1 = 244
 m2 = 42 
 sigma_P = 0.00069714285714285711
 sigma_L = 0.00012

 # Assignment rule definition
 !F S_V_Ratio = Mem_Area/Vcyt
 !F Mconc = (M)/M_init
 !F Lconc = (L)/L_init
 !F Pconc = (P)/P_init
 
 # Assignment rule initialisations
 M_init = 199693.0
 L_init = 102004
 P_init = 5303
 Mconc = 1.0
 Lconc = 1.0
 Pconc = 1.0

 # Species initialisations
 N@Vout = 3.07e5
 Pmem@Mem_Area = 37.38415
 Lmem@Mem_Area = 8291.2350678770199
 M@Vcyt = 199693.0
 L@Vcyt = 102004
 P@Vcyt = 5303
 
 # Parameter initialisations
 k1 = 0.00089709
 k2 = 0.000182027
 k3 = 1.7539e-010
 k4 = 5.0072346e-005
 k5 = 0.000574507164
 
 """
 Simulate this model to 200 for maximum happiness and
 watch the surface to volume ratio and scaled concentrations.
 """
 
This example illustrates almost all the new features included 
in the PySCeS MDL. Although it may be slightly more complicated 
than the basic model described above it is still, by our 
definition, human readable. 
