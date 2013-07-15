C
C*    Group  Linear Solver subroutines (Code DECCON/SOLCON)
C
      SUBROUTINE DECCON(A,NROW,NCOL,MCON,M,N,IRANKC,IRANK,COND,
     *                  D,PIVOT,KRED,AH,V,IERR)
C*    Begin Prologue DECCON
      INTEGER IRANKC,IRANK,MCON
      INTEGER M,N,NROW,NCOL,KRED
      INTEGER PIVOT(NCOL)
      DOUBLE PRECISION COND
      DOUBLE PRECISION A(NROW,NCOL),AH(NCOL,NCOL)
      DOUBLE PRECISION D(NCOL),V(NCOL)
      INTEGER IERR
C     ------------------------------------------------------------
C
C*  Title
C
C*    Deccon - Constrained Least Squares QR-Decomposition
C
C*  Written by        P. Deuflhard, U. Nowak, L. Weimann 
C*  Purpose           Solution of least squares problems, optionally
C                     with equality constraints.
C*  Method            Constrained Least Squares QR-Decomposition
C                     (see references below)
C*  Category          D9b1. -  Singular, overdetermined or
C                              underdetermined systems of linear 
C                              equations, generalized inverses. 
C                              Constrained Least Squares solution
C*  Keywords          Linear Least Square Problems, constrained, 
C                     QR-decomposition, pseudo inverse.
C*  Version           1.3
C*  Revision          December 1993
C*  Latest Change     August 2006
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad-Zuse-Zentrum fuer
C                     Informationstechnik Berlin (ZIB)
C                     Takustrasse 7, D-14195 Berlin-Dahlem
C                     phone : + 49/30/84185-0
C                     fax   : + 49/30/84185-125
C*  Contact           Lutz Weimann
C                     ZIB, Division Scientific Computing, 
C                          Department Numerical Analysis and Modelling
C                     phone : + 49/30/84185-185
C                     fax   : + 49/30/84185-107
C                     e-mail: weimann@zib.de
C
C*    References:
C     ===========
C
C       /1/ P.Deuflhard, V.Apostolescu:
C           An underrelaxed Gauss-Newton method for equality
C           constrained nonlinear least squares problems.
C           Lecture Notes Control Inform. Sci. vol. 7, p.
C           22-32 (1978)
C       /2/ P.Deuflhard, W.Sautter:
C           On rank-deficient pseudoinverses.
C           J. Lin. Alg. Appl. vol. 29, p. 91-111 (1980)
C    
C*    Related Programs:     SOLCON
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time. 
C    In any case you should not deliver this code without a special 
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty 
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from acquisition or application of this code.
C
C* Software status 
C    This code is under care of ZIB and belongs to ZIB software class 1.
C
C     ------------------------------------------------------------
C
C*    Summary:
C     ========
C     Constrained QR-decomposition of (M,N)-system  with
C     computation of pseudoinverse in case of rank-defeciency .
C     First MCON rows belong to equality constraints.
C
C     ------------------------------------------------------------
C
C*    Parameters list description (* marks inout parameters)
C     ======================================================
C
C*    Input parameters
C     ================
C
C       A(NROW,NCOL) Dble   Array holding the (M,N)-Matrix to be 
C                           decomposed
C       NROW         Int    Declared number of rows of array A
C       NCOL         Int    Declared number of columns of array A and 
C                           rows and columns of array AH
C       MCON         Int    Number of equality constraints (MCON.LE.N)
C                           Internally reduced if equality constraints
C                           are linearly dependent
C       M            Int    Current number of rows of matrix A
C       N            Int    Current number of columns of matrix A
C     * IRANKC       Int    Prescribed maximum pseudo-rank of 
C                           constrained part of matrix A (IRANKC.LE.MCON)
C     * IRANK        Int    Prescribed maximum pseudo-rank of matrix A
C                           (IRANK.LE.N)
C     * COND         Dble   Permitted upper bound for the subcondition
C                           of the least squares part of A, .i.e.
C                           DABS(D(IRANKC+1)/D(IRANK))
C       KRED         Int    Type of operation
C                           >=0  Householder triangularization
C                                (build up pseudo-inverse,if IRANK.LT.N)
C                           < 0  Reduction of pseudo-rank of matrix A, 
C                                skipping Householder triangularization,
C                                 build-up new pseudo-inverse
C
C*    Output parameters
C     =================
C
C       A(NROW,NCOL)  Dble   Array holding the (M,N)-output consisting
C                            of the transformed matrix in the upper 
C                            right triangle and the performed House-
C                            holder transf. in the lower left triangle.
C     * IRANKC        Int    New pseudo-rank of constrained part of
C                            matrix A, determined so that
C                            DABS(D(1)/D(IRANKC))<1/EPMACH
C     * IRANK         Int    New pseudo-rank of matrix A, determined
C                            so that DABS(D(IRANKC+1)/D(IRANK)) < COND
C       D(IRANK)      Dble   Diagonal elements of upper triangular matr.
C       PIVOT(N)      Int    Index vector storing permutation of columns
C                            due to pivoting
C     * COND          Dble   The sub-condition number belonging to the
C                            least squares part of A.
C                            (in case of rank reduction:
C                             sub-condition number which led to
C                             rank reduction)
C                            COND=0 indicates COND=infinity
C       AH(NCOL,NCOL) Dble   In case of rank-defect used to compute the
C                            pseudo-inverse (currently used will be an
C                            (N,N)-part of this array)
C       V(N)          Dble   V(1) holds on output the sub-condition
C                            number belonging to the constrained part
C                            of A.
C       IERR          Int    Error indicator:
C                            = 0 : DECCON computations are successfull.
C                            =-2 : Numerically negative diagonal element
C                                  encountered during computation of
C                                  pseudo inverse - due to extremely bad
C                                  conditioned Matrix A. DECCON is
C                                  unable to continue rank-reduction.
C
C*    Workspace parameters
C     ====================
C
C       V(N)         Dble   Workspace array
C
C*    Subroutines called: ZIBCONST
C
C*    Machine constants used
C     ======================
C
C     EPMACH = relative machine precision
      DOUBLE PRECISION EPMACH, SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS,DSQRT
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION REDUCE
      PARAMETER (REDUCE=0.05D0)
      INTEGER L1
      DOUBLE PRECISION S1
      INTEGER I,II,IRANKH,IRK1,I1,J,JD,JJ,K,K1,LEVEL,MH, IDATA
      DOUBLE PRECISION DD,D1MACH,H,HMAX,S,SH,T
C*    Begin
C     --------------------------------------------------------------
C     1 Initialization
      CALL ZIBCONST(EPMACH,SMALL)
      IF(IRANK.GT.N) IRANK = N
      IF(IRANK.GT.M) IRANK = M
C     --------------------------------------------------------------
C     1.1 Special case M=1 and N=1
      IF(M.EQ.1.AND.N.EQ.1)THEN
        PIVOT(1)=1
        D(1)=A(1,1)
        COND = ONE
        RETURN
      ENDIF
      IF(KRED.GE.0) THEN
C       ------------------------------------------------------------
C       1.1 Initialize pivot-array
        DO 11 J=1,N
          PIVOT(J)=J
11      CONTINUE
C       ------------------------------------------------------------
C       2. Constrained Householder triangularization    
        JD = 1
        IRANC1 = IRANKC + 1
        MH = MCON
        IRANKH = IRANKC
        IDATA = 0
        IF(MH.EQ.0) THEN
          IRANKH = IRANK
          MH = M
          IDATA = 1
        ENDIF
        IRK1 = IRANK
        DO 2  K=1,IRK1
2000      LEVEL = 1
          IF(K.NE.N) THEN
             K1 = K+1
C            DO (Until)
20           CONTINUE
             IF(JD.NE.0) THEN
               DO 201 J=K,N
                 S = ZERO
                 DO 2011 L1=K,MH
                    S = S+A(L1,J)**2
2011             CONTINUE
                 D(J)=S
201           CONTINUE
             ENDIF
C            ------------------------------------------------------
C            2.1 Column pivoting
             S1 = D(K)
             JJ = K
             DO 21 L1=K,N
               IF(D(L1).GT.S1) THEN
                 S1=D(L1)
                 JJ = L1
               ENDIF
21           CONTINUE
             H = D(JJ)
             IF(JD.EQ.1) HMAX = H/(DMAX1(1.0D1,COND*REDUCE))
             JD = 0
             IF(H.LT.HMAX) JD = 1
             IF(.NOT.(H.GE.HMAX)) GOTO 20
C            UNTIL ( expression - negated above)
C
             IF(JJ.NE.K) THEN
C              ------------------------------------------------------
C               2.2 Column interchange
                I = PIVOT(K)
                PIVOT(K)=PIVOT(JJ)
                PIVOT(JJ)=I
                D(JJ)=D(K)
                DO 221 L1=1,M
                   S1=A(L1,JJ)
                   A(L1,JJ)=A(L1,K)
                   A(L1,K)=S1
221             CONTINUE
             ENDIF
C          endif for k.ne.n case
           ENDIF
           H = ZERO
           DO 222 L1=K,MH
              H = H+A(L1,K)**2
222        CONTINUE
           T = DSQRT(H)
C          ----------------------------------------------------------
C          2.3.0 A-priori test on pseudo-rank
           IF  ( K.EQ.1 .OR. K.EQ.IRANC1 )  DD = T/COND
           IF(T.LE.DD .OR. K.GT.IRANKH) THEN
C            ------------------------------------------------------
C             2.3.1 Rank reduction
              IRANKH = K-1
              IF  (MH.NE.MCON .OR. IDATA.EQ.1)  THEN
                 IRANK = IRANKH
                 IF (IRANKC.EQ.IRANK) THEN
                    LEVEL = 4
                 ELSE
                    LEVEL = 3
                 ENDIF
              ELSE
                 IRANKC = IRANKH
                 IF  (IRANKC.NE.MCON) THEN
                    MH = M
                    IRANKH = IRANK
                    JD = 1
                    IDATA = 1
                    GOTO 2000
                 ELSE
                   STOP 'INTERNAL ERROR OF DECCON'
                 ENDIF
              ENDIF
           ENDIF
C
           IF (LEVEL.EQ.1) THEN
C             ------------------------------------------------------
C             2.4 Householder step
              S = A(K,K)
              T = -DSIGN(T,S)
              D(K)=T
C             By updating a(k,k) at this stage the 241 and 242 loop
C             must not be modified for l1=k.
              A(K,K)=S-T
              IF(K.NE.N) THEN
                 T = ONE/(H-S*T)
                 DO 24 J=K1,N
                    S = ZERO
                    DO 241 L1=K,MH
                          S = S+A(L1,K)*A(L1,J)
241                 CONTINUE
                    S = S*T
                    S1 = -S
                    IF(S.NE.0.D0) THEN
C                      Update the sub columns
                       DO 242 L1=K,M
                          A(L1,J) = A(L1,J)+A(L1,K)*S1
242                    CONTINUE
                    ENDIF
C                   Update sub column norms
                    D(J) = D(J)-A(K,J)**2
24               CONTINUE
                 IF(K.EQ.IRANKC) THEN
                    MH = M
                    JD = 1
                    IRANKH = IRANK
                  ENDIF
                  IF (K.EQ.IRK1) LEVEL = 3
              ELSE
                 LEVEL = 4
              ENDIF
C endif Householder step
           ENDIF
C        Exit Do 2 If ... 
           IF(LEVEL.GT.1) GOTO  2999
2        CONTINUE
C        ENDDO
2999     CONTINUE
      ELSE
         K = -1
         LEVEL = 3
      ENDIF
C     --------------------------------------------------------------
C     3 Rank-deficient pseudo-inverse
      IF(LEVEL.EQ.3) THEN
         IRK1 = IRANK+1
         DO 3 J=IRK1,N
            DO 31 II=1,IRANK
               I = IRK1-II
               S = A(I,J)
               IF(II.NE.1) THEN
                  SH = ZERO
                  DO 311 L1=I1,IRANK
                     SH=SH+A(I,L1)*V(L1)
311               CONTINUE
                  S = S-SH
               ENDIF
               I1 = I
               V(I)=S/D(I)
               AH(I,J)=V(I)
31          CONTINUE
            DO 32 I=IRK1,J
              S = ZERO
              DO 321 L1=1,I-1
                 S = S+AH(L1,I)*V(L1)
321           CONTINUE
              IF(I.NE.J)THEN
                 V(I)=-S/D(I)
                 AH(I,J)=-V(I)
              ENDIF
32         CONTINUE
           IF (S.GT.-ONE) THEN
              D(J)=DSQRT(S+ONE)
           ELSE 
              IERR=-2
              GOTO 999
           ENDIF
3        CONTINUE
      ENDIF
C    --------------------------------------------------------------
C     9 Exit
      IF (IRANKC.NE.0) THEN
        SH = D(IRANKC)
        IF(SH.NE.ZERO) SH = DABS(D(1)/SH)
      ELSE
        SH=ZERO
      ENDIF
      V(1) = SH
      IF (K.EQ.IRANK) T = D(IRANK)
      IF (IRANKC+1.LE.IRANK .AND. T.NE.ZERO) THEN
        S = DABS(D(IRANKC+1)/T)
      ELSE
        S = ZERO
      ENDIF
      COND = S
      IERR=0
999   CONTINUE
      RETURN
      END

      SUBROUTINE SOLCON(A,NROW,NCOL,MCON,M,N,X,B,IRANKC,IRANK,
     *                  D,PIVOT,KRED,AH,V)
C*    Begin Prologue SOLCON
      DOUBLE PRECISION A(NROW,NCOL),AH(NCOL,NCOL)
      DOUBLE PRECISION X(NCOL),B(NROW),D(NCOL),V(NCOL)
      INTEGER NROW,NCOL,MCON,M,N,IRANKC,IRANK,KRED
      INTEGER PIVOT(NCOL)
C     ------------------------------------------------------------
C    
C*    Summary
C     =======
C
C     Best constrained linear least squares solution of (M,N)-
C     system . First MCON rows comprise MCON equality constraints.
C     To be used in connection with subroutine DECCON
C     References:       See DECCON
C     Related Programs: DECCON
C    
C*    Parameters:
C     ===========
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       A(M,N), NROW, NCOL, M, N, MCON, IRANKC, IRANK,
C       D(N), PIVOT(N), AH(N,N), KRED
C                           See input- respective output-parameters
C                           description of subroutine DECCON
C     * B(M)         Dble   Right-hand side of linear system, if
C                           KRED.GE.0
C                           Right-hand side of upper linear system,
C                           if KRED.LT.0
C
C*    Output parameters
C     =================
C
C       X(N)         Dble   Best LSQ-solution of linear system
C       B(M)         Dble   Right-hand of upper trigular system
C                           (transformed right-hand side of linear
C                            system)
C
C*    Workspace parameters
C     ====================
C
C       V(N)         Dble   Workspace array
C
C     ------------------------------------------------------------
C*    End Prologue
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1,L2
      INTEGER I,II,I1,IRANC1,IRK1,J,JJ,J1,MH
      DOUBLE PRECISION S,SH
C*    Begin
C     ------------------------------------------------------------
C     1 Solution for pseudo-rank zero
      IF(IRANK.EQ.0)THEN
        DO 11 L1=1,N
          X(L1)=ZERO
11      CONTINUE
        RETURN
      ENDIF
      IF  ((IRANK.LE.IRANKC).AND.(IRANK.NE.N)) THEN
        IRANC1 = IRANKC + 1
        DO 12 L1=IRANC1,N
          V(L1) = ZERO
12      CONTINUE
      ENDIF
      IF(KRED.GE.0.AND.(M.NE.1.OR.N.NE.1))THEN
C       ----------------------------------------------------------
C       2 Constrained householder transformations of right-hand side
        MH = MCON
        IF(IRANKC.EQ.0) MH = M
        DO 21 J=1,IRANK
          S = ZERO
          DO 211 L1=J,MH
            S = S+A(L1,J)*B(L1)
211       CONTINUE
          S = S/(D(J)*A(J,J))
          DO 212 L1=J,M
            B(L1)=B(L1)+A(L1,J)*S
212       CONTINUE
          IF(J.EQ.IRANKC) MH = M
21      CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     3 Solution of upper triangular system
      IRK1 = IRANK+1
      DO 31 II=1,IRANK
        I = IRK1-II
        I1 = I+1
        S = B(I)
        IF(II.NE.1)THEN 
          SH = ZERO
          DO 311 L1=I1,IRANK 
            SH=SH+A(I,L1)*V(L1)
311       CONTINUE
          S = S-SH
        ENDIF
        V(I)=S/D(I)
31    CONTINUE
      IF((IRANK.NE.N).AND.(IRANK.NE.IRANKC)) THEN
C       ----------------------------------------------------------
C       3.2 Computation of the best constrained least squares-
C           solution
        DO 321 J=IRK1,N
          S = ZERO
          DO 3211 L1=1,J-1
            S = S+AH(L1,J)*V(L1)
3211      CONTINUE
          V(J)=-S/D(J)
321     CONTINUE
        DO 322 JJ=1,N
          J = N-JJ+1
          S = ZERO
          IF(JJ.NE.1) THEN
            DO 3221 L1=J1,N
              S=S+AH(J,L1)*V(L1)
3221        CONTINUE
          ENDIF
          IF(JJ.NE.1.AND.J.LE.IRANK) THEN
            V(J)=V(J)-S
          ELSE
            J1 = J
            V(J)=-(S+V(J))/D(J)
          ENDIF
322     CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     4 Back-permutation of solution components
      DO 4 L1=1,N
        L2 = PIVOT(L1)
        X(L2) = V(L1)
4     CONTINUE
      RETURN
      END
C*    End package
