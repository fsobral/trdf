      SUBROUTINE EASYTRDF(N,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF,EVALC,
     +     EVALJAC,EVALHC,F,FEAS,FCNT)

      use trdf, only : trdfsub

      IMPLICIT NONE

C     This subroutine calls the true TRDF subroutine using default
C     configuration parameters

C     SCALAR ARGUMENTS
      integer M,N,FCNT
      double precision F,FEAS

C     ARRAY ARGUMENTS
      double precision  X(N),XL(N),XU(N)
      logical CCODED(2),EQUATN(M),LINEAR(M)

C     LOCAL SCALARS
      logical OUTPUT
      integer MAXFCNT,NPT
      double precision RBEG,REND,XEPS

C     EXTERNAL SUBROUTINES
      external EVALF,EVALC,EVALJAC,EVALHC

C     NUMBER OF INTERPOLATION POITNS. 

      if ( N .le. 2 ) then
         NPT = 2 * N + 1
      else
         NPT = 2 * N + 3
      end if

C     SETS DEFAULT PARAMETERS
      
      MAXFCNT = 100000
      
      RBEG = 1.0D-1
      REND = 1.0D-4
      XEPS = 1.0D-8

      OUTPUT = .true.

C     CALLS THE TRUE SUBROUTINE

      CALL TRDFSUB(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF,EVALC,
     +             EVALJAC,EVALHC,MAXFCNT,RBEG,REND,XEPS,OUTPUT,
     +             F,FEAS,FCNT)

      END

C     ******************************************************************
C     ******************************************************************

      SUBROUTINE FULLTRDF(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF,
     +                    EVALC,EVALJAC,EVALHC,MAXFCNT,RBEG,REND,XEPS,
     +                    OUTPUT,F,FEAS,FCNT)     

      use trdf, only : trdfsub

      IMPLICIT NONE

C     This subroutine is just an interface between f77 and the f90
C     implementation of TRDF

C     SCALAR ARGUMENTS
      logical output
      integer m,maxfcnt,N,NPT,FCNT
      double precision F,FEAS,RBEG,REND,XEPS

C     ARRAY ARGUMENTS
      DOUBLE PRECISION  X(N),XL(N),XU(N)
      logical ccoded(2),equatn(m),linear(m)

C     EXTERNAL SUBROUTINES
      external EVALF,EVALC,EVALJAC,EVALHC


      CALL TRDFSUB(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF,EVALC,
     +             EVALJAC,EVALHC,MAXFCNT,RBEG,REND,XEPS,OUTPUT,F,
     +             FEAS,FCNT)

      END
