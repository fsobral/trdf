      PROGRAM PRINCIPAL
      
      IMPLICIT NONE

#include "tr_params.par"

C     LOCAL SCALARS

      integer I,M,NN,NPT
      double precision RELDIFF

C     LOCAL ARRAYS

      LOGICAL ccoded(2),equatn(MMAX),linear(MMAX)
      double precision XX(NMAX),L(NMAX),U(NMAX)

C     USER-DEFINED SUBROUTINES

      external calobjf,calcon,caljac,calhc

C     COMMON SCALARS

      integer N,NILI,NINL,NELI,NENL,NEX,NTP
      logical INDEX1(MMAX),INDEX2(MMAX),LXL(NMAX),LXU(NMAX)
      double precision G(MMAX),GG(NMAX * MMAX),X(NMAX),XL(NMAX),XU(NMAX)
      logical LEX
      double precision FEX

C     COMMON ARRAYS

      double precision XEX(NMAX * NMAX)

C     COMMON BLOCKS

      common/L1/N,NILI,NINL,NELI,NENL
      common/L2/X
      common/L3/G
      common/L5/GG
      common/L8/NTP
      common/L9/INDEX1
      common/L10/INDEX2
      common/L11/LXL
      common/L12/LXU
      common/L13/XL
      common/L14/XU
      COMMON/L20/LEX,NEX,FEX,XEX

C     LOCAL SCALARS

      character OPTM
      integer FCNT,MAXFCNT
      double precision C,F,FEAS,RBEG,REND,XEPS

C     WRITE(*,*) 'Number of the problem: '
      READ(*,*)NTP

      CALL CONV(1)

C     NUMBER OF VARIABLES
      nn = N

C     NUMBER OF INTERPOLATION POITNS. 

      if ( nn .le. 2 ) then
         NPT = 2 * nn + 1
      else
         NPT = 2 * nn + 3
      end if

C     INITIAL POINT AND BOX CONSTRAINTS.

      do i = 1,nn
         xx(i) = X(i)
      end do

      do i = 1,nn
         if ( LXL(i) ) then
            l(i) = XL(i)
         else
            l(i) = - 1.0D+20
         end if
      end do

      do i = 1,nn
         if ( LXU(i) ) then
            u(i) = XU(i)
         else
            u(i) = 1.0D+20
         end if
      end do

C     NUMBER OF CONSTRAINTS
      M = NILI + NINL + NELI + NENL

C     CONSTRAINTS

      do i = 1,NILI
         linear(i) = .true.
         equatn(i) = .false.
      end do

      do i = NILI + 1,NILI + NINL
         linear(i) = .false.
         equatn(i) = .false.
      end do

      do i = NILI + NINL + 1,NILI + NINL + NELI
         linear(i) = .true.
         equatn(i) = .true.
      end do

      do i = NILI + NINL + NELI + 1,M
         linear(i) = .false.
         equatn(i) = .true.
      end do

      do i = 1,M
         INDEX1(i) = .false.
         INDEX2(i) = .false.
      end do

C     CODED SUBROUTINES FOR CONSTRAINTS' DERIVATIVES

      ccoded(1) = .true.
      ccoded(2) = .false.

C     MAXIMUM NUMBER OF FUNCTION EVALUATIONS
      MAXFCNT = 100000

C     Some HS problems do not have derivatives of the constraints
      if ( NTP .eq. 348 .or. NTP .eq. 332 .or. NTP .eq. 365 .or. 
     +     NTP .eq. 362 .or. NTP .eq. 363 .or. NTP .eq. 364 .or.
     +     NTP .eq. 365 .or. NTP .eq. 366 .or. NTP .eq. 369 .or.
     +     NTP .eq. 390 .or. NTP .eq. 393 ) then
         ccoded(1) = .false.
      end if

C     CALLS THE ALGORITHM

      open(75,FILE='runhs.out')
      write(75,0020) NTP,N,NILI + NINL,NELI + NENL,1.0D20,1.0D20,
     +     1.0D20,-1
      close(75)

      MAXFCNT = 100000
      
      RBEG = 1.0D-1
      REND = 1.0D-4
      XEPS = 1.0D-8

      CALL FULLTRDF(NN,NPT,XX,L,U,M,EQUATN,LINEAR,CCODED,CALOBJF,CALCON,
     +              CALJAC,CALHC,MAXFCNT,RBEG,REND,XEPS,F,FEAS,FCNT)

      reldiff = (f - FEX) / max(1.0D0,abs(f),abs(FEX))

      optm = ' '
      if ( reldiff .le. 1.0D-01 .and. FEAS .le. XEPS ) then
         optm = '*'
      end if

      open(75,FILE='runhs.out')
      write(75,0020) NTP,N,NILI + NINL,NELI + NENL,FEX,F,FEAS,FCNT
      write(*,0021) NTP,N,NILI + NINL,NELI + NENL,FEX,F,FEAS,FCNT,
     +     optm
      close(75)

!     NON-EXECUTABLE STATEMENTS

 0020 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E15.8,1X,E15.8,1X,E15.8,1X,I15)
 0021 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E15.8,1X,E15.8,1X,E15.8,1X,I15,
     +     1X,A1)

      END PROGRAM PRINCIPAL

C     ******************************************************************
C     ******************************************************************

      subroutine calobjf(n,xx,f)

      implicit none

#include "tr_params.par"

!     SCALAR ARGUMENTS
      integer n
      double precision f

!     ARRAY ARGUMENTS
      double precision xx(n)

!     COMMON SCALARS
      double precision FX

!     COMMON ARRAYS
      double precision X(NMAX)

!     COMMON BLOCKS
      common/L2/X
      common/L6/FX

!     LOCAL SCALARS
      integer i

      do i = 1,n
         X(i) = xx(i)
      end do

      CALL CONV(2)

      f = FX

      end subroutine calobjf

C     ******************************************************************
C     ******************************************************************

      subroutine calcon(nn,xx,ind,c)

      implicit none

#include "tr_params.par"

!     SCALAR ARGUMENTS
      integer ind,nn
      double precision c

!     ARRAY ARGUMENTS
      double precision xx(nn)

C     COMMON ARRAYS
      double precision G(MMAX),X(NMAX)
      logical INDEX1(MMAX)

C     COMMON BLOCKS
      common/L2/X
      common/L3/G
      common/L9/INDEX1

C     LOCAL SCALARS
      integer i

      do i = 1,nn
         X(i) = xx(i)
      end do
      
      INDEX1(ind) = .true.
      
      call conv(4)

      INDEX1(ind) = .false.
      
      c = - G(ind)

      end subroutine calcon

C     ******************************************************************
C     ******************************************************************

      subroutine caljac(nn,xx,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

      implicit none

#include "tr_params.par"

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,ind,jcnnz,lim,nn

C     ARRAY ARGUMENTS
      integer jcvar(lim)
      double precision xx(nn),jcval(lim)

C     COMMON SCALARS
      integer N,NILI,NINL,NELI,NENL

C     COMMON ARRAYS
      double precision GG(NMAX * MMAX),X(NMAX)
      logical INDEX2(MMAX)

C     COMMON BLOCKS
      common/L1/N,NILI,NINL,NELI,NENL
      common/L2/X
      common/L5/GG
      common/L10/INDEX2

C     LOCAL SCALARS
      integer i,m
      
      flag = 0

      lmem = .false.

      do i = 1,nn
         X(i) = xx(i)
      end do

      INDEX2(ind) = .true.
      
      CALL CONV(5)

      INDEX2(ind) = .false.
      
      jcnnz = nn

      m = NILI + NINL + NELI + NENL

      do i = 1,nn
         jcvar(i) = i
         jcval(i) = - GG((i - 1) * m + ind)
      end do 

      end subroutine caljac

C     ******************************************************************
C     ******************************************************************

      subroutine calhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,hcnnz,ind,lim,n

C     ARRAY ARGUMENTS
      integer hccol(lim),hcrow(lim)
      double precision hcval(lim),x(n)

      flag = -1

      lmem = .false.

      end
