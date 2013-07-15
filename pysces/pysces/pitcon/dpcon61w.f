C dpcon61w.f 29 February 2004

C PySCeS - Python Simulator for Cellular Systems 
C            (http://pysces.sourceforge.net)
C Copyright (C) B.G. Olivier, J.M. Rohwer, J.-H.S. Hofmeyr
C Stellenbosch, 2004-2008.
C Triple-J Group for Molecular Cell Physiology
C Stellenbosch University, South Africa
C Author: Brett G. Olivier
C
C PySCeS is Open Source Software distributed under
C the GNU GENERAL PUBLIC LICENSE (see docs/GPL)

CFILE: dpcon61w.f
      SUBROUTINE PITCON1(DF,FPAR,FX,IERROR,IPAR,IWORK,LIW,NVAR,RWORK,
     *LRW,XR,IMTH)
      EXTERNAL  DF
      EXTERNAL  FX
      EXTERNAL  DENSLV
      EXTERNAL  BANSLV
      INTEGER   LIW
cf2py integer optional,check(len(iwork)>=liw),depend(iwork) :: liw=len(iwork)
      INTEGER   LRW
cf2py integer optional,check(len(rwork)>=lrw),depend(rwork) :: lrw=len(rwork)
      INTEGER   NVAR
cf2py integer optional,check(len(xr)>=nvar),depend(xr) :: nvar=len(xr)
      DOUBLE PRECISION FPAR(*)
cf2py double precision dimension(*)),intent(in) :: fpar 
      INTEGER   IPAR(*)
cf2py integer dimension(*),intent(in) :: ipar 
      INTEGER   IWORK(LIW)
cf2py integer dimension(liw),intent(in,out,copy) :: iwork
      INTEGER   IERROR
cf2py integer intent(out) :: ierror
      DOUBLE PRECISION RWORK(LRW)
cf2py double precision dimension(lrw),intent(in,out,copy) :: rwork
      DOUBLE PRECISION XR(NVAR)
cf2py double precision dimension(nvar),intent(in,out,copy) :: xr
cf2py DOUBLE PRECISION X(NVAR), FVEC(NVAR), FJAC(NVAR,NVAR)
cf2py double precision,intent(out) :: fvec
cf2py double precision,intent(out) :: fjac
cf2py INTEGER IMTH = 1
      INTEGER IMTH
cf2py CALL FX(NVAR,FPAR,IPAR,X,FVEC,IERROR)
cf2py CALL DF(NVAR,FPAR,IPAR,X,FJAC,IERROR)
      IF (IMTH.eq.1) THEN
         CALL PITCON(DF,FPAR,FX,IERROR,IPAR,IWORK,LIW,NVAR,RWORK,
     *LRW,XR,DENSLV)
      ELSE
         CALL PITCON(DF,FPAR,FX,IERROR,IPAR,IWORK,LIW,NVAR,RWORK,
     *LRW,XR,BANSLV)
      ENDIF
      END
