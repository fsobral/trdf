      PROGRAM PRINCIPAL

#include "tr_params.par"

C     LOCAL ARRAYS
      LOGICAL ccoded(2),equatn(MMAX),linear(MMAX)
      double precision X(NMAX),XL(NMAX),XU(NMAX)

C     USER-DEFINED SUBROUTINES
      external calobjf,calcon,caljac,calhc

C     NUMBER OF VARIABLES
      N= 3

C     NUMBER OF CONSTRAINTS
      M = 2

C     TYPE OF CONSTRAINTS
      do i = 1,M
         equatn(i) = .false.
      end do
      
      do i = i,m
         linear(i) = .true.
      end do

C     INITIAL POINT AND BOX CONSTRAINTS.
      DO 6 I=1,3
         X(I)=10.D0
         XL(I)= 0D0
 6       XU(I)= 42D0 

C     CODED SUBROUTINES FOR CONSTRAINTS' DERIVATIVES
      ccoded(1) = .true.
      ccoded(2) = .true.

C     CALLS THE ALGORITHM

      CALL EASYTRDF(N,X,XL,XU,M,EQUATN,LINEAR,CCODED,CALOBJF,CALCON,
     +              CALJAC,CALHC,F,FEAS,FCNT)

      END PROGRAM PRINCIPAL

C     ******************************************************************
C     ******************************************************************

      subroutine calobjf(n,x,f,flag)

!     SCALAR ARGUMENTS
      integer flag,n
      double precision f

!     ARRAY ARGUMENTS
      double precision x(n)

!     Problem HS37 from Hock-Schittkowski collection

      flag = 0

      f = - x(1) * x(2) * x(3)

      end subroutine calobjf

C     ******************************************************************
C     ******************************************************************

      subroutine calcon(n,x,ind,c,flag)

!     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

!     ARRAY ARGUMENTS
      double precision x(n)

!     Problem HS37 from Hock-Schittkowski collection

      flag = 0

      if ( ind .eq. 1 ) then
         c = - 72.D0 + x(1) + 2.D0 * x(2) + 2.D0 * x(3)
         return
      elseif ( ind .eq. 2 ) then
         c = - x(1) - 2.D0 * x(2) - 2.D0 * x(3)
         return
      end if

      flag = -1

      end subroutine calcon

C     ******************************************************************
C     ******************************************************************

      subroutine caljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,ind,jcnnz,lim,n

C     ARRAY ARGUMENTS
      integer jcvar(lim)
      double precision x(n),jcval(lim)

      lmem = .false.

      flag = 0

!     Problem HS37 from Hock-Schittkowski collection

      if ( ind .eq. 1 ) then
         jcnnz    = 3
         jcvar(1) = 1
         jcval(1) = 1.D0  
         jcvar(2) = 2
         jcval(2) = 2.D0  
         jcvar(3) = 3
         jcval(3) = 2.D0  
         return
      else if ( ind .eq. 2 ) then
         jcnnz    = 3
         jcvar(1) = 1
         jcval(1) = -1.D0
         jcvar(2) = 2
         jcval(2) = -2.D0     
         jcvar(3) = 3
         jcval(3) = -2.D0     
         return
      else
         
         flag = -1

      end if
  
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

!     Problem HS37 from Hock-Schittkowski collection

      flag = 0

      lmem = .false.

!     The problem has linear constraints

      hcnnz = 0

      end
