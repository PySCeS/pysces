NLEQ2 package - release 2.3 at January 3, 1992

Description of the package:
---------------------------

The NLEQ2 package consists of the files listed and briefly described
below:

filename    contents of file

nleq2.f     NLEQ2 standard user interface and internal subroutines.

wnorm.f     Norm computation subroutine used for termination criteria.

main_nleq2.f    An example main program illustrating the usage of NLEQ2.

nleq2e.f    An easy-to-use driver subroutine to NLEQ2 for nonlinear
            problems with n <= 50 equations. The calling parameter list
            consists of 4 items and, additionally, only the nonlinear
            problem function (named FCN) must be supplied.

main_nleq2_easy.f    An example main program illustrating the usage of NLEQ2E.

linalg_nleq2.f Linear algebra routines (decoon, solcon).

zibmon.f    Monitor.

zibsec.f    Time routine ZIBSEC.

zibconst.f   Machine constants ZIBCONST.

makefile    An input file to the UNIX make utility containing informa-
            tion how to build the executable programs nleq2 and
            nleq2-easy.

readme      This information.

nleq2_out.nrm   Example monitor output of program run of nleq2.

nleq2_dat.nrm   Example data output of program run of nleq2.

nleq2_easy.nrm  Example monitor output (to FORTRAN unit 6) of program run
            of nleq2-easy.


The two main programs solve the same sequence of test problems 
(the Chebyquad problem with dimensions n=2,3,...,9) in two
slightly different ways. To build the executable programs, you need to
compile and to link:

main_nleq2.f, nleq2.f, linalg_nleq2.f, zibmon.f, zibsec.f,
  zibconst.f   (to build program nleq2),
main_nleq2_easy.f, nleq2e.f, nleq2.f, linalg_nleq2.f, zibmon.f,
  zibsec.f, zibconst.f   (to build program nleq2-easy).

Under UNIX with the make-utility, for example simply type in:
make             -   to build program nleq2    - or
make nleq2-easy  -   to build program nleq2-easy.

Note that a few adaptations may be necessary to utilize the code for
your computer at hand. Please examine the two subroutines ZIBCONST and
ZIBSEC. Probably, you need to adapt ZIBCONST to correspond to the arithmetic 
of your computer. If you want to obtain time measurements from the time
monitor, you must supply a "time stamp", measured in seconds, through the 
argument of the subroutine zibsec.
If you replace version 2.2 with the new one and if you examine the
output values from RWK (e.g. CONV, SUMX, DLEVF), please notice that
the positions of some RWK output values have been changed.

What is new within the current release?
---------------------------------------

A lot of improvements have been included within release 2.3 of NLEQ2
since the previous release 2.2. They are shortly summarized below:

- The easy to use interface NLEQ2E - see above
- The time monitor, which can be used to measure the performance
  of the package and of the users problem subroutines -
  see IOPT(19), IOPT(20) and subroutine ZIBSEC
- An optionally selectable bounded damping strategy which can be helpful
  for solving extremely nonlinear problems -
  see IOPT(38) (also IOPT(31)) and RWK(20)
- If the user problem function FCN returns IERR=1, indicating that
  F(X1) cannot be evaluated for the given X1, a damping factor reduction
  is done and F is tried to be evaluated again with a different argument
  X2. The degree of reduction of the damping factor may be also
  controlled by the user - see description of FCN
- Certain characteristic internal parameters of the numerical approxi-
  mation subroutines may be altered now easily by the user -
  see RWK(26), RWK(27), RWK(28), RWK(29)
- The convergence rate of the Newton iteration is monitored and based
  in this, the iteration will be terminated if it significantly slows
  down when approaching the solution - due to roundoff errors. (A
  corresponding warning message will be printed and a warning code
  returned through parameter IERR). Also, in connection with this 
  monitor feature, a warning message and code will be given, if 
  the iteration terminated because the usual criterion has been 
  fulfilled, but the corresponding error estimate may not be appro-
  priate due to lack of theoretical assumptions in the actual context.
  see IOPT(39) and error messages
- The standard rank-1 updates procedure, operating directly on the
  Jacobian, has been replaced by an (equivalent) iterative rank-1
  updates algorithm which uses exclusively vectors for computations -
  see IOPT(32), IWK(36).
- Double precision storage will be saved against the previous version
  if rank-1 updates are inhibited - up to n(n-1) elements.
- The computation of the norm which is used in connection with the
  convergence termination criterion, is now separated in the function
  subroutine WNORM and may be easily exchanged.
- NLEQ2 is now a real extension of NLEQ1 applied with full mode Jacobian
  - in the sense that the damping factor calculations produce the same
  sequence of damping factors as NLEQ1 - until damping factor reduction
  would lead to error termination of NLEQ1 - the case when rank 
  reduction is activated by NLEQ2. (but practically iteration sequences
  of NLEQ1 and NLEQ2 may differ before this condition meets because
  of using different linear solver procedures)

Further minor changes:

- The default value for IOPT(3) (JACGEN) is now 2 (Jacobian by numerical
  approximation  w i t h o u t  feedback device).
- RWK positions, which have been changed:
  name     new   old
  CONV     17    24
  SUMX     18    26
  DLEVF    19    27
  SIGMA2   24    28



 
