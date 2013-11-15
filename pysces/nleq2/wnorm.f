      DOUBLE PRECISION FUNCTION WNORM(N,Z,XW)
      INTEGER N
      DOUBLE PRECISION Z(N), XW(N)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     E N O R M : Return the norm to be used in exit (termination)
C                 criteria
C
C*    Input parameters
C     ================
C
C     N         Int Number of equations/unknowns
C     Z(N)     Dble  The vector, of which the norm is to be computed
C     XW(N)    Dble  The scaling values of Z(N)
C
C*    Output
C     ======
C
C     WNORM(N,Z,XW)  Dble  The mean square root norm of Z(N) subject
C                          to the scaling values in XW(N):
C                          = Sqrt( Sum(1,...N)((Z(I)/XW(I))**2) / N )
C
C     ------------------------------------------------------------
C*    End Prologue
      INTEGER I
      DOUBLE PRECISION S
C*    Begin
      S = 0.0D0
      DO 10 I=1,N
        S = S + ( Z(I)/XW(I) ) ** 2
10    CONTINUE
      WNORM = DSQRT( S / DBLE(FLOAT(N)) )
C     End of function WNORM
      RETURN
      END
